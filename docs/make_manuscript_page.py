#!/usr/bin/env python3
"""make_manuscript_page.py MANUSCRIPT.md FIGURE.png OUT.html
Renders the manuscript draft (markdown) as one self-contained HTML page for the artifact viewer:
[AUTHOR ...] and [NOT RUN ...] markers are highlighted, the chaos-to-order sketch is embedded after its paragraph."""
import base64, re, sys, markdown
md, fig, out = sys.argv[1:4]
def split_paragraphs(text):          # the draft writes one paragraph per line: put a blank line between plain adjacent lines
    L, out = text.split('\n'), []
    plain = lambda l: l.strip() and not re.match(r'\s*([-*] |\d+\. |\||#)', l) or re.match(r'\s*\*[A-Za-z]', l) is not None
    for i, l in enumerate(L):
        out.append(l)
        if i + 1 < len(L) and l.strip() and L[i + 1].strip() and plain(l) and plain(L[i + 1]) and not l.startswith('|'): out.append('')
    return '\n'.join(out)
body = markdown.markdown(split_paragraphs(open(md, encoding='utf-8').read()), extensions=['tables', 'sane_lists'])
body = re.sub(r'\[(AUTHOR[^\]]*|NOT RUN[^\]]*)\]', lambda m: '<mark>[%s]</mark>' % m.group(1), body)
img = base64.b64encode(open(fig, 'rb').read()).decode()
figure = ('<figure><img alt="Sketch of Figure 2: 390 real B2 copies in random order, in SubFam order, and the 20 chunk consensuses" '
          'src="data:image/png;base64,%s"><figcaption>Figure 2 (sketch). 390 real B2 copies (stratified random sample): input order, SubFam order, 20 chunk consensuses on one alignment. '
          'Left strip: Vassetzky 2021 best-hit label. Red: differs from the column majority; white: gap.</figcaption></figure>' % img)
body = re.sub(r'(<p><em>From chaos to order.*?</p>)', lambda m: m.group(1) + figure, body, count=1, flags=re.S)
title = 'SubFam manuscript draft'
css = open(__file__.replace('make_manuscript_page.py', 'manuscript_page.css'), encoding='utf-8').read()
open(out, 'w', encoding='utf-8').write('<title>%s</title>\n<style>\n%s\n</style>\n<main>\n%s\n</main>\n' % (title, css, body))
