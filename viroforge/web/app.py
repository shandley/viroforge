#!/usr/bin/env python3
"""
Flask Application for ViroForge Web Interface

Author: ViroForge Development Team
Date: 2025-11-10
"""

import os
from pathlib import Path
from flask import Flask, render_template, request, jsonify, send_file, redirect, url_for
from werkzeug.utils import secure_filename
import json
import subprocess
import tempfile

# Import ViroForge utilities
from viroforge.database.utils import get_all_collections, get_collection_details
from viroforge.cli.presets import list_presets as get_presets, load_preset
from viroforge.cli.report import load_dataset_metadata, load_composition_file
from viroforge.cli.batch import load_batch_config, expand_parameter_sweep

app = Flask(__name__)
app.config['SECRET_KEY'] = 'viroforge-web-interface-secret-key'
app.config['MAX_CONTENT_LENGTH'] = 16 * 1024 * 1024  # 16MB max upload

# Store running generations (in production, use Redis or similar)
running_generations = {}


def create_app():
    """Create and configure the Flask application."""
    return app


@app.route('/')
def index():
    """Home page."""
    return render_template('index.html')


@app.route('/collections')
def collections():
    """Browse collections."""
    try:
        all_collections = get_all_collections()
        return render_template('collections.html', collections=all_collections)
    except Exception as e:
        return render_template('error.html', error=str(e)), 500


@app.route('/api/collections')
def api_collections():
    """API: Get all collections."""
    try:
        collections = get_all_collections()
        return jsonify(collections)
    except Exception as e:
        return jsonify({'error': str(e)}), 500


@app.route('/api/collection/<int:collection_id>')
def api_collection_details(collection_id):
    """API: Get collection details."""
    try:
        details = get_collection_details(collection_id)
        if details is None:
            return jsonify({'error': 'Collection not found'}), 404
        return jsonify(details)
    except Exception as e:
        return jsonify({'error': str(e)}), 500


@app.route('/generate')
def generate_page():
    """Dataset generation page."""
    try:
        presets = get_presets()
        collections = get_all_collections()
        return render_template('generate.html', presets=presets, collections=collections)
    except Exception as e:
        return render_template('error.html', error=str(e)), 500


@app.route('/api/presets')
def api_presets():
    """API: Get all presets."""
    try:
        presets = get_presets()
        return jsonify(presets)
    except Exception as e:
        return jsonify({'error': str(e)}), 500


@app.route('/api/preset/<preset_name>')
def api_preset_details(preset_name):
    """API: Get preset details."""
    try:
        preset = load_preset(preset_name)
        if preset is None:
            return jsonify({'error': 'Preset not found'}), 404
        return jsonify(preset)
    except Exception as e:
        return jsonify({'error': str(e)}), 500


@app.route('/api/generate', methods=['POST'])
def api_generate():
    """API: Start dataset generation."""
    try:
        data = request.get_json()

        # Build viroforge generate command
        cmd = ['viroforge', 'generate']

        if data.get('preset'):
            cmd.extend(['--preset', data['preset']])
        else:
            # Direct parameters
            if data.get('collection_id'):
                cmd.extend(['--collection-id', str(data['collection_id'])])
            if data.get('platform'):
                cmd.extend(['--platform', data['platform']])
            if data.get('coverage'):
                cmd.extend(['--coverage', str(data['coverage'])])
            if data.get('depth'):
                cmd.extend(['--depth', str(data['depth'])])

        # Override parameters
        if data.get('output'):
            cmd.extend(['--output', data['output']])
        if data.get('seed'):
            cmd.extend(['--seed', str(data['seed'])])
        if data.get('vlp_protocol'):
            cmd.extend(['--vlp-protocol', data['vlp_protocol']])
        if data.get('no_vlp'):
            cmd.append('--no-vlp')
        if data.get('contamination_level'):
            cmd.extend(['--contamination-level', data['contamination_level']])

        # Start generation in background
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        # Store process
        generation_id = f"gen_{process.pid}"
        running_generations[generation_id] = {
            'process': process,
            'command': ' '.join(cmd),
            'status': 'running'
        }

        return jsonify({
            'generation_id': generation_id,
            'command': ' '.join(cmd),
            'status': 'started'
        })

    except Exception as e:
        return jsonify({'error': str(e)}), 500


@app.route('/api/generation/<generation_id>/status')
def api_generation_status(generation_id):
    """API: Check generation status."""
    if generation_id not in running_generations:
        return jsonify({'error': 'Generation not found'}), 404

    gen = running_generations[generation_id]
    process = gen['process']

    if process.poll() is None:
        # Still running
        return jsonify({
            'status': 'running',
            'command': gen['command']
        })
    else:
        # Completed
        return_code = process.returncode
        stdout, stderr = process.communicate()

        status = 'success' if return_code == 0 else 'failed'
        gen['status'] = status

        return jsonify({
            'status': status,
            'return_code': return_code,
            'stdout': stdout,
            'stderr': stderr,
            'command': gen['command']
        })


@app.route('/batch')
def batch_page():
    """Batch generation page."""
    return render_template('batch.html')


@app.route('/api/batch/validate', methods=['POST'])
def api_batch_validate():
    """API: Validate batch YAML configuration."""
    try:
        yaml_content = request.get_json().get('yaml')

        # Save to temp file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            f.write(yaml_content)
            temp_path = f.name

        try:
            # Try to load configuration
            config = load_batch_config(temp_path)

            # Count datasets
            datasets = config.get('datasets', [])
            if 'parameter_sweep' in config:
                sweep_datasets = expand_parameter_sweep(config['parameter_sweep'])
                datasets.extend(sweep_datasets)

            return jsonify({
                'valid': True,
                'batch_name': config.get('batch_name'),
                'total_datasets': len(datasets),
                'datasets': [d.get('name', f'dataset_{i}') for i, d in enumerate(datasets, 1)]
            })

        finally:
            os.unlink(temp_path)

    except Exception as e:
        return jsonify({
            'valid': False,
            'error': str(e)
        })


@app.route('/api/batch/generate', methods=['POST'])
def api_batch_generate():
    """API: Start batch generation."""
    try:
        yaml_content = request.get_json().get('yaml')
        parallel = request.get_json().get('parallel', 1)

        # Save to temp file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            f.write(yaml_content)
            temp_path = f.name

        # Build command
        cmd = ['viroforge', 'batch', temp_path]
        if parallel > 1:
            cmd.extend(['--parallel', str(parallel)])

        # Start batch generation
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        # Store process
        generation_id = f"batch_{process.pid}"
        running_generations[generation_id] = {
            'process': process,
            'command': ' '.join(cmd),
            'status': 'running',
            'temp_file': temp_path
        }

        return jsonify({
            'generation_id': generation_id,
            'command': ' '.join(cmd),
            'status': 'started'
        })

    except Exception as e:
        return jsonify({'error': str(e)}), 500


@app.route('/reports')
def reports_page():
    """Dataset reports page."""
    return render_template('reports.html')


@app.route('/api/report', methods=['POST'])
def api_report():
    """API: Get dataset report."""
    try:
        dataset_path = Path(request.get_json().get('dataset_path'))

        if not dataset_path.exists():
            return jsonify({'error': 'Dataset not found'}), 404

        # Load metadata
        metadata = load_dataset_metadata(dataset_path)
        if not metadata:
            return jsonify({'error': 'Could not load metadata'}), 404

        # Load composition
        composition = load_composition_file(dataset_path)

        # Prepare response
        report = {
            'dataset_name': dataset_path.name,
            'metadata': metadata,
            'composition': composition if composition else []
        }

        return jsonify(report)

    except Exception as e:
        return jsonify({'error': str(e)}), 500


@app.route('/compare')
def compare_page():
    """Dataset comparison page."""
    return render_template('compare.html')


@app.route('/api/compare', methods=['POST'])
def api_compare():
    """API: Compare datasets."""
    try:
        dataset_paths = request.get_json().get('datasets', [])

        if len(dataset_paths) < 2:
            return jsonify({'error': 'Need at least 2 datasets to compare'}), 400

        # Load metadata for all datasets
        datasets = []
        for path_str in dataset_paths:
            dataset_path = Path(path_str)
            if not dataset_path.exists():
                continue

            metadata = load_dataset_metadata(dataset_path)
            if metadata:
                metadata['_path'] = str(dataset_path)
                metadata['_name'] = dataset_path.name
                datasets.append(metadata)

        if len(datasets) < 2:
            return jsonify({'error': 'Could not load metadata for datasets'}), 404

        # Analyze datasets
        collection_ids = set(d.get('collection', {}).get('id') for d in datasets)
        seeds = set(d.get('generation_info', {}).get('random_seed') for d in datasets if d.get('generation_info'))
        platforms = set(d.get('platform', {}).get('name') for d in datasets if d.get('platform'))

        # Determine comparison type
        comparison_type = 'unknown'
        recommendations = []

        if len(collection_ids) == 1 and len(seeds) == 1 and len(platforms) > 1:
            comparison_type = 'technology'
            recommendations.append({
                'type': 'success',
                'title': 'Suitable for technology/platform comparison',
                'details': [
                    'Same collection and seed ensures identical genome composition',
                    'Multiple platforms enable direct performance comparison'
                ]
            })

        # Check for hybrid assembly
        has_short = any(d.get('platform', {}).get('name') in ['novaseq', 'miseq', 'hiseq'] for d in datasets)
        has_long = any(d.get('platform', {}).get('name') in ['pacbio-hifi', 'nanopore'] for d in datasets)

        if has_short and has_long and len(collection_ids) == 1 and len(seeds) == 1:
            comparison_type = 'hybrid'
            recommendations.append({
                'type': 'success',
                'title': 'Suitable for hybrid assembly!',
                'details': [
                    'Short + long reads with matched compositions',
                    'Try: Unicycler, SPAdes hybrid mode, or MaSuRCA'
                ]
            })

        elif len(collection_ids) > 1:
            comparison_type = 'multi_collection'
            recommendations.append({
                'type': 'info',
                'title': 'Multi-collection comparison',
                'details': [
                    'Compare virome characteristics across different environments',
                    'Note: Different genomes, so assembly metrics not directly comparable'
                ]
            })

        elif len(platforms) == 1:
            comparison_type = 'same_platform'
            recommendations.append({
                'type': 'info',
                'title': 'Same-platform comparison',
                'details': [
                    'Useful for testing different parameters (coverage, VLP protocol, etc.)'
                ]
            })

        return jsonify({
            'datasets': datasets,
            'collection_ids': list(collection_ids),
            'seeds': list(seeds),
            'platforms': list(platforms),
            'comparison_type': comparison_type,
            'recommendations': recommendations
        })

    except Exception as e:
        return jsonify({'error': str(e)}), 500


@app.route('/about')
def about_page():
    """About page."""
    return render_template('about.html')


if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5000)
