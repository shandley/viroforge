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
import math
import subprocess
import tempfile

# Import ViroForge utilities
from viroforge.cli.db_utils import load_all_collections as get_all_collections, load_collection_details as get_collection_details
from viroforge.cli.preset_loader import list_all_presets as get_presets, load_preset
from viroforge.cli.report import load_dataset_metadata, load_composition_file
from viroforge.cli.batch import load_batch_config, expand_parameter_sweep

app = Flask(__name__)
app.config['SECRET_KEY'] = 'viroforge-web-interface-secret-key'
app.config['MAX_CONTENT_LENGTH'] = 16 * 1024 * 1024  # 16MB max upload


def _is_db_missing(exc):
    """Check if an exception is caused by an unpopulated database."""
    msg = str(exc).lower()
    return 'no such table' in msg or 'unable to open database' in msg

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
        if _is_db_missing(e):
            return render_template('setup_required.html'), 503
        return render_template('error.html', error=str(e)), 500


@app.route('/api/collections')
def api_collections():
    """API: Get all collections."""
    try:
        collections = get_all_collections()
        return jsonify(collections)
    except Exception as e:
        if _is_db_missing(e):
            return jsonify({'error': 'Database not set up. Run: viroforge setup-db'}), 503
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
        if _is_db_missing(e):
            return jsonify({'error': 'Database not set up. Run: viroforge setup-db'}), 503
        return jsonify({'error': str(e)}), 500


@app.route('/generate')
def generate_page():
    """Dataset generation page."""
    try:
        presets_by_category = get_presets()
        # Flatten grouped presets into a single list for the template
        presets = [p for category_presets in presets_by_category.values() for p in category_presets]
        collections = get_all_collections()
        return render_template('generate.html', presets=presets, collections=collections)
    except Exception as e:
        if _is_db_missing(e):
            return render_template('setup_required.html'), 503
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
        if data.get('collection_id'):
            cmd.extend(['--collection-id', str(data['collection_id'])])
        if data.get('platform'):
            cmd.extend(['--platform', data['platform']])
        if data.get('coverage'):
            cmd.extend(['--coverage', str(data['coverage'])])
        if data.get('depth'):
            cmd.extend(['--depth', str(data['depth'])])
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
        if data.get('dark_matter_fraction') is not None and data.get('dark_matter_fraction') != '':
            cmd.extend(['--dark-matter-fraction', str(data['dark_matter_fraction'])])

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
            'status': 'running',
            'output_dir': data.get('output', ''),
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

        report_url = None
        if status == 'success' and gen.get('output_dir'):
            output_path = Path(gen['output_dir'])
            report_file = output_path / f"{output_path.name}_report.html"
            if report_file.is_file():
                report_url = f"/report-file/{output_path.name}"
                gen['report_path'] = str(report_file)

        return jsonify({
            'status': status,
            'return_code': return_code,
            'report_url': report_url,
            'stdout': stdout,
            'stderr': stderr,
            'command': gen['command']
        })


@app.route('/report-file/<dataset_name>')
def serve_report(dataset_name):
    """Serve an auto-generated HTML report for a dataset."""
    for gen in running_generations.values():
        if gen.get('report_path') and Path(gen['report_path']).name.startswith(dataset_name):
            return send_file(gen['report_path'], mimetype='text/html')
    return jsonify({'error': 'Report not found'}), 404


@app.route('/api/report/generate', methods=['POST'])
def api_generate_report():
    """API: Generate HTML report for an existing dataset."""
    try:
        dataset_path = Path(request.get_json().get('dataset_path', '')).resolve()
        if not dataset_path.is_dir():
            return jsonify({'error': 'Dataset directory not found'}), 404
        from viroforge.cli.auto_report import generate_auto_report
        report_path = generate_auto_report(dataset_path)
        if report_path is None:
            return jsonify({'error': 'No metadata found in dataset'}), 404
        return send_file(str(report_path.resolve()), mimetype='text/html')
    except Exception as e:
        return jsonify({'error': str(e)}), 500


@app.route('/batch')
def batch_page():
    """Batch generation page."""
    try:
        collections = get_all_collections()
    except Exception:
        collections = []
    return render_template('batch.html', collections=collections)


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

        # Load composition and sanitize NaN values (not valid JSON)
        composition = load_composition_file(dataset_path) or []
        for entry in composition:
            for key, val in entry.items():
                if isinstance(val, float) and math.isnan(val):
                    entry[key] = None

        # Prepare response
        report = {
            'dataset_name': dataset_path.name,
            'metadata': metadata,
            'composition': composition
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
        # v1.1: platform at configuration.platform; v1.0: platform at platform.name
        def _get_platform(d):
            return d.get('platform', {}).get('name') or d.get('configuration', {}).get('platform')
        platforms = set(_get_platform(d) for d in datasets if _get_platform(d))

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
        has_short = any(_get_platform(d) in ['novaseq', 'miseq', 'hiseq'] for d in datasets)
        has_long = any(_get_platform(d) in ['nanopore'] for d in datasets)

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


@app.route('/benchmark')
def benchmark_page():
    """Benchmark page."""
    return render_template('benchmark.html')


@app.route('/api/benchmark', methods=['POST'])
def api_benchmark():
    """API: Start a benchmark run."""
    try:
        data = request.get_json()
        module = data.get('module')
        if module not in ('qc', 'assembly', 'taxonomy'):
            return jsonify({'error': 'Invalid module'}), 400

        cmd = ['viroforge', 'benchmark', module]

        if module == 'qc':
            raw = data.get('raw_reads', '')
            cleaned = data.get('cleaned_reads', '')
            if not raw or not cleaned:
                return jsonify({'error': 'raw_reads and cleaned_reads are required'}), 400
            cmd.extend(['--raw-reads', raw])
            cmd.extend(['--cleaned-reads', cleaned])
            if data.get('keep_remove'):
                for pair in data['keep_remove'].split():
                    cmd.extend(['--keep-remove', pair])

        elif module == 'assembly':
            contigs = data.get('contigs', '')
            genomes = data.get('genomes', '')
            if not contigs or not genomes:
                return jsonify({'error': 'contigs and genomes are required'}), 400
            cmd.extend(['--contigs', contigs])
            cmd.extend(['--genomes', genomes])
            if data.get('ground_truth'):
                cmd.extend(['--ground-truth', data['ground_truth']])

        elif module == 'taxonomy':
            pipeline_output = data.get('pipeline_output', '')
            ground_truth = data.get('ground_truth', '')
            if not ground_truth:
                return jsonify({'error': 'ground_truth is required'}), 400
            if pipeline_output:
                cmd.extend(['--pipeline-output', pipeline_output])
            cmd.extend(['--ground-truth', ground_truth])
            fmt = data.get('format', 'auto')
            if fmt and fmt != 'auto':
                cmd.extend(['--format', fmt])
            mode = data.get('mode', 'read-based')
            if mode != 'read-based':
                cmd.extend(['--mode', mode])
            if mode == 'contig-based':
                if data.get('contig_taxonomy'):
                    cmd.extend(['--contig-taxonomy', data['contig_taxonomy']])
                if data.get('contigs'):
                    cmd.extend(['--contigs', data['contigs']])
                if data.get('genomes'):
                    cmd.extend(['--genomes', data['genomes']])
            if data.get('taxdump_dir'):
                cmd.extend(['--taxdump-dir', data['taxdump_dir']])

        output_json = tempfile.NamedTemporaryFile(
            mode='w', suffix='.json', delete=False, prefix='bench_'
        )
        output_json.close()
        markdown_file = output_json.name.replace('.json', '.md')
        cmd.extend(['--output', output_json.name])
        cmd.extend(['--markdown', markdown_file])

        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        benchmark_id = f"bench_{process.pid}"
        running_generations[benchmark_id] = {
            'process': process,
            'command': ' '.join(cmd),
            'status': 'running',
            'module': module,
            'output_json': output_json.name,
            'markdown_file': markdown_file,
        }

        return jsonify({
            'benchmark_id': benchmark_id,
            'command': ' '.join(cmd),
            'status': 'started'
        })

    except Exception as e:
        return jsonify({'error': str(e)}), 500


@app.route('/api/benchmark/<benchmark_id>/status')
def api_benchmark_status(benchmark_id):
    """API: Check benchmark status."""
    if benchmark_id not in running_generations:
        return jsonify({'error': 'Benchmark not found'}), 404

    bench = running_generations[benchmark_id]
    process = bench['process']

    if process.poll() is None:
        return jsonify({'status': 'running', 'command': bench['command']})

    return_code = process.returncode
    stdout, stderr = process.communicate()
    status = 'success' if return_code == 0 else 'failed'
    bench['status'] = status

    result = {
        'status': status,
        'return_code': return_code,
        'stdout': stdout,
        'stderr': stderr,
        'command': bench['command'],
    }

    if status == 'success':
        try:
            with open(bench['output_json']) as f:
                result['metrics'] = json.load(f)
        except Exception:
            pass
        try:
            with open(bench['markdown_file']) as f:
                result['markdown'] = f.read()
        except Exception:
            pass
        if result.get('metrics'):
            try:
                from viroforge.benchmarking.html_report import generate_benchmark_report
                html_path = bench['output_json'].replace('.json', '_report.html')
                generate_benchmark_report(result['metrics'], bench['module'], html_path)
                bench['report_path'] = html_path
                result['report_url'] = f'/benchmark-report/{benchmark_id}'
            except Exception:
                pass

    return jsonify(result)


@app.route('/benchmark-report/<benchmark_id>')
def serve_benchmark_report(benchmark_id):
    """Serve the HTML benchmark report."""
    if benchmark_id not in running_generations:
        return jsonify({'error': 'Benchmark not found'}), 404
    bench = running_generations[benchmark_id]
    report_path = bench.get('report_path')
    if not report_path or not Path(report_path).is_file():
        return jsonify({'error': 'Report not found'}), 404
    return send_file(report_path, mimetype='text/html')


@app.route('/about')
def about_page():
    """About page."""
    return render_template('about.html')


if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5000)
