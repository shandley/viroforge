#!/usr/bin/env python3
"""
Web Command

Launch ViroForge web interface.

Author: ViroForge Development Team
Date: 2025-11-10
"""

import sys
import webbrowser
import time
from threading import Timer

try:
    from rich.console import Console
except ImportError:
    print("Error: 'rich' library required. Install with: pip install rich", file=sys.stderr)
    sys.exit(1)

console = Console()


def run_web(args):
    """Launch web interface."""
    try:
        from viroforge.web.app import app
    except ImportError as e:
        console.print("[red]Error: Flask not installed[/red]")
        console.print()
        console.print("The web interface requires Flask. Install with:")
        console.print("  [cyan]pip install flask[/cyan]")
        console.print()
        console.print(f"Error details: {e}")
        return 1

    host = args.host
    port = args.port
    debug = args.debug

    console.print()
    console.print("═" * 80, style="bold cyan")
    console.print(" ViroForge Web Interface", style="bold cyan")
    console.print("═" * 80, style="bold cyan")
    console.print()
    console.print(f"Starting server on [cyan]http://{host}:{port}[/cyan]")
    console.print()
    console.print("Features:")
    console.print("  • Browse virome collections")
    console.print("  • Generate datasets with presets")
    console.print("  • Build batch configurations")
    console.print("  • View dataset reports")
    console.print("  • Compare multiple datasets")
    console.print()
    console.print("[dim]Press Ctrl+C to stop the server[/dim]")
    console.print()

    # Open browser after short delay
    if not args.no_browser:
        def open_browser():
            webbrowser.open(f'http://{host}:{port}')

        Timer(1.5, open_browser).start()

    try:
        app.run(host=host, port=port, debug=debug)
    except KeyboardInterrupt:
        console.print()
        console.print("[yellow]Server stopped[/yellow]")
        return 0
    except Exception as e:
        console.print(f"[red]Error starting server: {e}[/red]")
        return 1

    return 0
