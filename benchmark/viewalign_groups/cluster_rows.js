#!/usr/bin/env node
'use strict';
/* cluster_rows.js MSAVIEWER_DIR IN.aln.fasta OUT.json [--trim] [--min-occupancy F] [--min-perfect N] [--min-occ N] [--max-iter N] [--min-size N]
 * Runs ViewAlign's position-pattern subfamily clustering (cluster.js, class SINEClusterer, "Cluster Now") on an alignment, exactly the
 * code of the checkout MSAVIEWER_DIR, loaded the way tests/clustering/*.test.js load it. Defaults are the UI defaults used in those tests
 * (minSize 3, minPerfect 5, maxIterations 10, quality 80/70/60 %, sizes 11/20, minOccurrences 3). Writes the clusters (member row ids,
 * perfect = exclusive columns, imperfect, cloudy) and the unassigned rows as JSON.
 * --recurse D re-runs the clustering inside every cluster of at least 2*minSize rows, up to depth D (columns are kept, so positions stay the same):
 * each cluster then has 'children'. --min-perfect applies to every level. */
const fs = require('fs'), path = require('path'), vm = require('vm');
const a = process.argv.slice(2); const root = a[0], inp = a[1], out = a[2];
const opt = (f, d) => { const i = a.indexOf(f); return i >= 0 ? +a[i + 1] : d; };
const src = fs.readFileSync(path.join(root, 'cluster.js'), 'utf8');
const ctx = { console: { log() {}, warn() {}, error: console.error }, performance, setTimeout, MessageChannel, Map, Set, Array, Math, Object };
vm.createContext(ctx); vm.runInContext(src + '\nthis.SINEClusterer = SINEClusterer; this.getTrimBoundaries = getTrimBoundaries;', ctx);
const seqs = []; let cur = null;
for (const l of fs.readFileSync(inp, 'utf8').split(/\r?\n/)) {
    if (l.startsWith('>')) { if (cur) seqs.push(cur); cur = { id: l.slice(1).trim().split(/\s+/)[0], seq: '' }; } else if (cur) cur.seq += l.trim().toUpperCase();
}
if (cur) seqs.push(cur);
const UI = { minSize: opt('--min-size', 3), minPerfect: opt('--min-perfect', 5), maxIterations: opt('--max-iter', 10), qualitySmall: 80, qualityMedium: 70, qualityLarge: 60,
             sizeSmallMedium: 11, sizeMediumLarge: 20, minOccurrences: opt('--min-occ', 3) };
if (opt('--min-occupancy', 0) > 0) {   // mask sparse columns (all rows set to '-'): a base seen in 3 rows of a column that is a gap in 100 rows is not 'exclusive' evidence
    const th = opt('--min-occupancy', 0), L = seqs[0].seq.length, keep = new Array(L).fill(false);
    for (let j = 0; j < L; j++) { let n = 0; for (const x of seqs) if (x.seq[j] !== '-' && x.seq[j] !== '.') n++; keep[j] = n / seqs.length >= th; }
    for (const x of seqs) x.seq = x.seq.split('').map((c, j) => keep[j] ? c : '-').join('');
}
if (a.includes('--trim')) { const t = ctx.getTrimBoundaries(seqs); UI.trimStart = t.leftTrimEnd; UI.trimEnd = t.rightTrimStart; }
const slim = c => ({ size: c.size, rows: c.sequences.map(s => s.id), perfect: (c.perfectFeatures || []).map(f => f.pos + ':' + f.char),
                     imperfect: (c.imperfectFeatures || []).map(f => f.pos + ':' + f.char), cloudy: (c.cloudyFeatures || []).map(f => f.pos + ':' + f.char) });
async function run(rows, depth) {
    const r = await new ctx.SINEClusterer(rows).clusterChunked(UI);
    const out = r.clusters.map(slim);
    if (depth > 0) for (const c of out) {
        if (c.size >= 2 * UI.minSize && c.size < rows.length) { const ids = new Set(c.rows); const sub = await run(rows.filter(x => ids.has(x.id)), depth - 1); if (sub.clusters.length > 1) c.children = sub.clusters; c.childUnassigned = sub.unassigned; }
    }
    return { clusters: out, unassigned: r.unassigned.map(s => s.id), summary: r.summary };
}
if (a.includes('--groups')) {          // --groups ROWS_GROUPS.tsv: characterise a GIVEN partition (row<TAB>group), no discovery: strict = exclusive columns, cloudy = enriched
    const g = {}; for (const l of fs.readFileSync(a[a.indexOf('--groups') + 1], 'utf8').split(/\r?\n/)) { const f = l.split('\t'); if (f.length >= 2) (g[f[1]] = g[f[1]] || []).push(f[0]); }
    const names = Object.keys(g).sort(), idx = new Map(seqs.map((x, i) => [x.id, i]));
    const ch = new ctx.SINEClusterer(seqs).characterizeGroups(names.map(n => g[n].map(id => idx.get(id)).filter(i => i !== undefined)), { minSize: UI.minSize });
    const res = {}; names.forEach((n, i) => res[n] = { rows: g[n].length, strict: ch[i].perfectFeatures.map(f => f.pos + ':' + f.char), cloudy: ch[i].cloudyFeatures.map(f => f.pos + ':' + f.char) });
    fs.writeFileSync(out, JSON.stringify(res, null, 1));
    for (const n of names) console.log('%s: %d rows, %d strict, %d cloudy', n, res[n].rows, res[n].strict.length, res[n].cloudy.length);
    process.exit(0);
}
(async () => {
    const r = await run(seqs, opt('--recurse', 0));
    fs.writeFileSync(out, JSON.stringify({ options: UI, nRows: seqs.length, alnLen: seqs[0].seq.length, clusters: r.clusters, unassigned: r.unassigned }, null, 1));
    console.log('%s: %d rows x %d columns -> %d top-level clusters, %d rows assigned, %d unassigned', path.basename(inp), seqs.length, seqs[0].seq.length, r.clusters.length, r.summary.nAssigned, r.unassigned.length);
    process.exit(0);   // the MessageChannel of cluster.js keeps node alive otherwise
})();
