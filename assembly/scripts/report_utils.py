#!/usr/bin/env python3
"""Write Markdown and escaped standalone HTML reports."""
import html
from pathlib import Path


def write_reports(markdown, html_path, text):
    for path in [Path(markdown), Path(html_path)]:
        path.parent.mkdir(parents=True, exist_ok=True)
    Path(markdown).write_text(text)
    # Self-contained, escaped HTML; no remote JavaScript or untrusted Markdown rendering.
    Path(html_path).write_text('<!doctype html><html lang="zh-CN"><meta charset="utf-8">'
                               '<title>Assembly report</title><style>body{max-width:1000px;'
                               'margin:32px auto;font:16px/1.6 sans-serif}pre{white-space:pre-wrap;'
                               'overflow-wrap:anywhere}</style><body><pre>'
                               + html.escape(text) + '</pre></body></html>\n')
