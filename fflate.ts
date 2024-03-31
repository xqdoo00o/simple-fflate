const u8 = Uint8Array, u16 = Uint16Array, i32 = Int32Array;

// fixed length extra bits
const fleb = new u8([0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0, /* unused */ 0, 0, /* impossible */ 0]);

// fixed distance extra bits
const fdeb = new u8([0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, /* unused */ 0, 0]);

// code length index map
const clim = new u8([16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15]);

// get base, reverse index map from extra bits
const freb = (eb: Uint8Array, start: number) => {
    const b = new u16(31);
    for (let i = 0; i < 31; ++i) {
        b[i] = start += 1 << eb[i - 1];
    }
    // numbers here are at max 18 bits
    const r = new i32(b[30]);
    for (let i = 1; i < 30; ++i) {
        for (let j = b[i]; j < b[i + 1]; ++j) {
            r[j] = ((j - b[i]) << 5) | i;
        }
    }
    return { b, r };
}

const { b: fl, r: revfl } = freb(fleb, 2);
// we can ignore the fact that the other numbers are wrong; they never happen anyway
fl[28] = 258, revfl[258] = 28;
const { b: fd, r: revfd } = freb(fdeb, 0);

// map of value to reverse (assuming 16 bits)
const rev = new u16(32768);
for (let i = 0; i < 32768; ++i) {
    // reverse table algorithm from SO
    let x = ((i & 0xAAAA) >> 1) | ((i & 0x5555) << 1);
    x = ((x & 0xCCCC) >> 2) | ((x & 0x3333) << 2);
    x = ((x & 0xF0F0) >> 4) | ((x & 0x0F0F) << 4);
    rev[i] = (((x & 0xFF00) >> 8) | ((x & 0x00FF) << 8)) >> 1;
}

// create huffman tree from u8 "map": index -> code length for code index
// mb (max bits) must be at most 15
// TODO: optimize/split up?
const hMap = ((cd: Uint8Array, mb: number, r: 0 | 1) => {
    const s = cd.length;
    // index
    let i = 0;
    // u16 "map": index -> # of codes with bit length = index
    const l = new u16(mb);
    // length of cd must be 288 (total # of codes)
    for (; i < s; ++i) {
        if (cd[i]) ++l[cd[i] - 1];
    }
    // u16 "map": index -> minimum code for bit length = index
    const le = new u16(mb);
    for (i = 1; i < mb; ++i) {
        le[i] = (le[i - 1] + l[i - 1]) << 1;
    }
    let co: Uint16Array;
    if (r) {
        // u16 "map": index -> number of actual bits, symbol for code
        co = new u16(1 << mb);
        // bits to remove for reverser
        const rvb = 15 - mb;
        for (i = 0; i < s; ++i) {
            // ignore 0 lengths
            if (cd[i]) {
                // num encoding both symbol and bits read
                const sv = (i << 4) | cd[i];
                // free bits
                const r = mb - cd[i];
                // start value
                let v = le[cd[i] - 1]++ << r;
                // m is end value
                for (const m = v | ((1 << r) - 1); v <= m; ++v) {
                    // every 16 bit value starting with the code yields the same result
                    co[rev[v] >> rvb] = sv;
                }
            }
        }
    } else {
        co = new u16(s);
        for (i = 0; i < s; ++i) {
            if (cd[i]) {
                co[i] = rev[le[cd[i] - 1]++] >> (15 - cd[i]);
            }
        }
    }
    return co;
});

// fixed length tree
const flt = new u8(288);
for (let i = 0; i < 144; ++i) flt[i] = 8;
for (let i = 144; i < 256; ++i) flt[i] = 9;
for (let i = 256; i < 280; ++i) flt[i] = 7;
for (let i = 280; i < 288; ++i) flt[i] = 8;
// fixed distance tree
const fdt = new u8(32);
for (let i = 0; i < 32; ++i) fdt[i] = 5;
// fixed length map
const flm = /*#__PURE__*/ hMap(flt, 9, 0), flrm = /*#__PURE__*/ hMap(flt, 9, 1);
// fixed distance map
const fdm = /*#__PURE__*/ hMap(fdt, 5, 0), fdrm = /*#__PURE__*/ hMap(fdt, 5, 1);

// find max of array
const max = (a: Uint8Array | number[]) => {
    let m = a[0];
    for (let i = 1; i < a.length; ++i) {
        if (a[i] > m) m = a[i];
    }
    return m;
};

// read d, starting at bit p and mask with m
const bits = (d: Uint8Array, p: number, m: number) => {
    const o = (p / 8) | 0;
    return ((d[o] | (d[o + 1] << 8)) >> (p & 7)) & m;
}

// read d, starting at bit p continuing for at least 16 bits
const bits16 = (d: Uint8Array, p: number) => {
    const o = (p / 8) | 0;
    return ((d[o] | (d[o + 1] << 8) | (d[o + 2] << 16)) >> (p & 7));
}

// get end of byte
const shft = (p: number) => ((p + 7) / 8) | 0;

// typed array slice - allows garbage collector to free original reference,
// while being more compatible than .slice
const slc = (v: Uint8Array, s: number, e?: number) => {
    if (s == null || s < 0) s = 0;
    if (e == null || e > v.length) e = v.length;
    // can't use .constructor in case user-supplied
    return new u8(v.subarray(s, e));
}

// inflate state
type InflateState = {
    // lmap
    l?: Uint16Array;
    // dmap
    d?: Uint16Array;
    // lbits
    m?: number;
    // dbits
    n?: number;
    // final
    f?: number;
    // pos
    p?: number;
    // byte
    b?: number;
    // lstchk
    i: number;
};

// error codes
const ec = [
    'unexpected EOF',
    'invalid block type',
    'invalid length/literal',
    'invalid distance',
    'stream finished',
    'no stream handler',
    , // determined by compression function
    'no callback',
    'invalid UTF-8 data',
    'extra field too long',
    'date not in range 1980-2099',
    'filename too long',
    'stream finishing',
    'invalid zip data'
    // determined by unknown compression method
];

/**
 * An error generated within this library
 */
interface FlateError extends Error {
    /**
     * The code associated with this error
     */
    code: number;
};

const err = (ind: number, msg?: string | 0, nt?: 1) => {
    const e: Partial<FlateError> = new Error(msg || ec[ind]);
    e.code = ind;
    if ((Error as any).captureStackTrace) (Error as any).captureStackTrace(e, err);
    if (!nt) throw e;
    return e as FlateError;
}

// expands raw DEFLATE data
const inflt = (dat: Uint8Array, st: InflateState, buf?: Uint8Array, dict?: Uint8Array) => {
    // source length       dict length
    const sl = dat.length, dl = dict ? dict.length : 0;
    if (!sl || st.f && !st.l) return buf || new u8(0);
    const noBuf = !buf;
    // have to estimate size
    const resize = noBuf || st.i != 2;
    // no state
    const noSt = st.i;
    // Assumes roughly 33% compression ratio average
    if (noBuf) buf = new u8(sl * 3);
    // ensure buffer can fit at least l elements
    const cbuf = (l: number) => {
        let bl = buf.length;
        // need to increase size to fit
        if (l > bl) {
            // Double or set to necessary, whichever is greater
            const nbuf = new u8(Math.max(bl * 2, l));
            nbuf.set(buf);
            buf = nbuf;
        }
    };
    //  last chunk         bitpos           bytes
    let final = st.f || 0, pos = st.p || 0, bt = st.b || 0, lm = st.l, dm = st.d, lbt = st.m, dbt = st.n;
    // total bits
    const tbts = sl * 8;
    do {
        if (!lm) {
            // BFINAL - this is only 1 when last chunk is next
            final = bits(dat, pos, 1);
            // type: 0 = no compression, 1 = fixed huffman, 2 = dynamic huffman
            const type = bits(dat, pos + 1, 3);
            pos += 3;
            if (!type) {
                // go to end of byte boundary
                const s = shft(pos) + 4, l = dat[s - 4] | (dat[s - 3] << 8), t = s + l;
                if (t > sl) {
                    if (noSt) err(0);
                    break;
                }
                // ensure size
                if (resize) cbuf(bt + l);
                // Copy over uncompressed data
                buf.set(dat.subarray(s, t), bt);
                // Get new bitpos, update byte count
                st.b = bt += l, st.p = pos = t * 8, st.f = final;
                continue;
            }
            else if (type == 1) lm = flrm, dm = fdrm, lbt = 9, dbt = 5;
            else if (type == 2) {
                //  literal                            lengths
                const hLit = bits(dat, pos, 31) + 257, hcLen = bits(dat, pos + 10, 15) + 4;
                const tl = hLit + bits(dat, pos + 5, 31) + 1;
                pos += 14;
                // length+distance tree
                const ldt = new u8(tl);
                // code length tree
                const clt = new u8(19);
                for (let i = 0; i < hcLen; ++i) {
                    // use index map to get real code
                    clt[clim[i]] = bits(dat, pos + i * 3, 7);
                }
                pos += hcLen * 3;
                // code lengths bits
                const clb = max(clt), clbmsk = (1 << clb) - 1;
                // code lengths map
                const clm = hMap(clt, clb, 1);
                for (let i = 0; i < tl;) {
                    const r = clm[bits(dat, pos, clbmsk)];
                    // bits read
                    pos += r & 15;
                    // symbol
                    const s = r >> 4;
                    // code length to copy
                    if (s < 16) {
                        ldt[i++] = s;
                    } else {
                        //  copy   count
                        let c = 0, n = 0;
                        if (s == 16) n = 3 + bits(dat, pos, 3), pos += 2, c = ldt[i - 1];
                        else if (s == 17) n = 3 + bits(dat, pos, 7), pos += 3;
                        else if (s == 18) n = 11 + bits(dat, pos, 127), pos += 7;
                        while (n--) ldt[i++] = c;
                    }
                }
                //    length tree                 distance tree
                const lt = ldt.subarray(0, hLit), dt = ldt.subarray(hLit);
                // max length bits
                lbt = max(lt)
                // max dist bits
                dbt = max(dt);
                lm = hMap(lt, lbt, 1);
                dm = hMap(dt, dbt, 1);
            } else err(1);
            if (pos > tbts) {
                if (noSt) err(0);
                break;
            }
        }
        // Make sure the buffer can hold this + the largest possible addition
        // Maximum chunk size (practically, theoretically infinite) is 2^17
        if (resize) cbuf(bt + 131072);
        const lms = (1 << lbt) - 1, dms = (1 << dbt) - 1;
        let lpos = pos;
        for (; ; lpos = pos) {
            // bits read, code
            const c = lm[bits16(dat, pos) & lms], sym = c >> 4;
            pos += c & 15;
            if (pos > tbts) {
                if (noSt) err(0);
                break;
            }
            if (!c) err(2);
            if (sym < 256) buf[bt++] = sym;
            else if (sym == 256) {
                lpos = pos, lm = null;
                break;
            } else {
                let add = sym - 254;
                // no extra bits needed if less
                if (sym > 264) {
                    // index
                    const i = sym - 257, b = fleb[i];
                    add = bits(dat, pos, (1 << b) - 1) + fl[i];
                    pos += b;
                }
                // dist
                const d = dm[bits16(dat, pos) & dms], dsym = d >> 4;
                if (!d) err(3);
                pos += d & 15;
                let dt = fd[dsym];
                if (dsym > 3) {
                    const b = fdeb[dsym];
                    dt += bits16(dat, pos) & (1 << b) - 1, pos += b;
                }
                if (pos > tbts) {
                    if (noSt) err(0);
                    break;
                }
                if (resize) cbuf(bt + 131072);
                const end = bt + add;
                if (bt < dt) {
                    const shift = dl - dt, dend = Math.min(dt, end);
                    if (shift + bt < 0) err(3);
                    for (; bt < dend; ++bt) buf[bt] = dict[shift + bt];
                }
                for (; bt < end; ++bt) buf[bt] = buf[bt - dt];
            }
        }
        st.l = lm, st.p = lpos, st.b = bt, st.f = final;
        if (lm) final = 1, st.m = lbt, st.d = dm, st.n = dbt;
    } while (!final)
    // don't reallocate for streams or user buffers
    return bt != buf.length && noBuf ? slc(buf, 0, bt) : buf.subarray(0, bt);
}

// starting at p, write the minimum number of bits that can hold v to d
const wbits = (d: Uint8Array, p: number, v: number) => {
    v <<= p & 7;
    const o = (p / 8) | 0;
    d[o] |= v;
    d[o + 1] |= v >> 8;
}

// starting at p, write the minimum number of bits (>8) that can hold v to d
const wbits16 = (d: Uint8Array, p: number, v: number) => {
    v <<= p & 7;
    const o = (p / 8) | 0;
    d[o] |= v;
    d[o + 1] |= v >> 8;
    d[o + 2] |= v >> 16;
}

type HuffNode = {
    // symbol
    s: number;
    // frequency
    f: number;
    // left child
    l?: HuffNode;
    // right child
    r?: HuffNode;
};

// creates code lengths from a frequency table
const hTree = (d: Uint16Array, mb: number) => {
    // Need extra info to make a tree
    const t: HuffNode[] = [];
    for (let i = 0; i < d.length; ++i) {
        if (d[i]) t.push({ s: i, f: d[i] });
    }
    const s = t.length;
    const t2 = t.slice();
    if (!s) return { t: et, l: 0 };
    if (s == 1) {
        const v = new u8(t[0].s + 1);
        v[t[0].s] = 1;
        return { t: v, l: 1 };
    }
    t.sort((a, b) => a.f - b.f);
    // after i2 reaches last ind, will be stopped
    // freq must be greater than largest possible number of symbols
    t.push({ s: -1, f: 25001 });
    let l = t[0], r = t[1], i0 = 0, i1 = 1, i2 = 2;
    t[0] = { s: -1, f: l.f + r.f, l, r };
    // efficient algorithm from UZIP.js
    // i0 is lookbehind, i2 is lookahead - after processing two low-freq
    // symbols that combined have high freq, will start processing i2 (high-freq,
    // non-composite) symbols instead
    // see https://reddit.com/r/photopea/comments/ikekht/uzipjs_questions/
    while (i1 != s - 1) {
        l = t[t[i0].f < t[i2].f ? i0++ : i2++];
        r = t[i0 != i1 && t[i0].f < t[i2].f ? i0++ : i2++];
        t[i1++] = { s: -1, f: l.f + r.f, l, r };
    }
    let maxSym = t2[0].s;
    for (let i = 1; i < s; ++i) {
        if (t2[i].s > maxSym) maxSym = t2[i].s;
    }
    // code lengths
    const tr = new u16(maxSym + 1);
    // max bits in tree
    let mbt = ln(t[i1 - 1], tr, 0);
    if (mbt > mb) {
        // more algorithms from UZIP.js
        // TODO: find out how this code works (debt)
        //  ind    debt
        let i = 0, dt = 0;
        //    left            cost
        const lft = mbt - mb, cst = 1 << lft;
        t2.sort((a, b) => tr[b.s] - tr[a.s] || a.f - b.f);
        for (; i < s; ++i) {
            const i2 = t2[i].s;
            if (tr[i2] > mb) {
                dt += cst - (1 << (mbt - tr[i2]));
                tr[i2] = mb;
            } else break;
        }
        dt >>= lft;
        while (dt > 0) {
            const i2 = t2[i].s;
            if (tr[i2] < mb) dt -= 1 << (mb - tr[i2]++ - 1);
            else ++i;
        }
        for (; i >= 0 && dt; --i) {
            const i2 = t2[i].s;
            if (tr[i2] == mb) {
                --tr[i2];
                ++dt;
            }
        }
        mbt = mb;
    }
    return { t: new u8(tr), l: mbt };
}
// get the max length and assign length codes
const ln = (n: HuffNode, l: Uint16Array, d: number): number => {
    return n.s == -1
        ? Math.max(ln(n.l, l, d + 1), ln(n.r, l, d + 1))
        : (l[n.s] = d);
}

// length codes generation
const lc = (c: Uint8Array) => {
    let s = c.length;
    // Note that the semicolon was intentional
    while (s && !c[--s]);
    const cl = new u16(++s);
    //  ind      num         streak
    let cli = 0, cln = c[0], cls = 1;
    const w = (v: number) => { cl[cli++] = v; }
    for (let i = 1; i <= s; ++i) {
        if (c[i] == cln && i != s)
            ++cls;
        else {
            if (!cln && cls > 2) {
                for (; cls > 138; cls -= 138) w(32754);
                if (cls > 2) {
                    w(cls > 10 ? ((cls - 11) << 5) | 28690 : ((cls - 3) << 5) | 12305);
                    cls = 0;
                }
            } else if (cls > 3) {
                w(cln), --cls;
                for (; cls > 6; cls -= 6) w(8304);
                if (cls > 2) w(((cls - 3) << 5) | 8208), cls = 0;
            }
            while (cls--) w(cln);
            cls = 1;
            cln = c[i];
        }
    }
    return { c: cl.subarray(0, cli), n: s };
}

// calculate the length of output from tree, code lengths
const clen = (cf: Uint16Array, cl: Uint8Array) => {
    let l = 0;
    for (let i = 0; i < cl.length; ++i) l += cf[i] * cl[i];
    return l;
}

// writes a fixed block
// returns the new bit pos
const wfblk = (out: Uint8Array, pos: number, dat: Uint8Array) => {
    // no need to write 00 as type: TypedArray defaults to 0
    const s = dat.length;
    const o = shft(pos + 2);
    out[o] = s & 255;
    out[o + 1] = s >> 8;
    out[o + 2] = out[o] ^ 255;
    out[o + 3] = out[o + 1] ^ 255;
    for (let i = 0; i < s; ++i) out[o + i + 4] = dat[i];
    return (o + 4 + s) * 8;
}

// writes a block
const wblk = (dat: Uint8Array, out: Uint8Array, final: number, syms: Int32Array, lf: Uint16Array, df: Uint16Array, eb: number, li: number, bs: number, bl: number, p: number) => {
    wbits(out, p++, final);
    ++lf[256];
    const { t: dlt, l: mlb } = hTree(lf, 15);
    const { t: ddt, l: mdb } = hTree(df, 15);
    const { c: lclt, n: nlc } = lc(dlt);
    const { c: lcdt, n: ndc } = lc(ddt);
    const lcfreq = new u16(19);
    for (let i = 0; i < lclt.length; ++i) ++lcfreq[lclt[i] & 31];
    for (let i = 0; i < lcdt.length; ++i) ++lcfreq[lcdt[i] & 31];
    const { t: lct, l: mlcb } = hTree(lcfreq, 7);
    let nlcc = 19;
    for (; nlcc > 4 && !lct[clim[nlcc - 1]]; --nlcc);
    const flen = (bl + 5) << 3;
    const ftlen = clen(lf, flt) + clen(df, fdt) + eb;
    const dtlen = clen(lf, dlt) + clen(df, ddt) + eb + 14 + 3 * nlcc + clen(lcfreq, lct) + 2 * lcfreq[16] + 3 * lcfreq[17] + 7 * lcfreq[18];
    if (bs >= 0 && flen <= ftlen && flen <= dtlen) return wfblk(out, p, dat.subarray(bs, bs + bl));
    let lm: Uint16Array, ll: Uint8Array, dm: Uint16Array, dl: Uint8Array;
    wbits(out, p, 1 + (dtlen < ftlen as unknown as number)), p += 2;
    if (dtlen < ftlen) {
        lm = hMap(dlt, mlb, 0), ll = dlt, dm = hMap(ddt, mdb, 0), dl = ddt;
        const llm = hMap(lct, mlcb, 0);
        wbits(out, p, nlc - 257);
        wbits(out, p + 5, ndc - 1);
        wbits(out, p + 10, nlcc - 4);
        p += 14;
        for (let i = 0; i < nlcc; ++i) wbits(out, p + 3 * i, lct[clim[i]]);
        p += 3 * nlcc;
        const lcts = [lclt, lcdt];
        for (let it = 0; it < 2; ++it) {
            const clct = lcts[it];
            for (let i = 0; i < clct.length; ++i) {
                const len = clct[i] & 31;
                wbits(out, p, llm[len]), p += lct[len];
                if (len > 15) wbits(out, p, (clct[i] >> 5) & 127), p += clct[i] >> 12;
            }
        }
    } else {
        lm = flm, ll = flt, dm = fdm, dl = fdt;
    }
    for (let i = 0; i < li; ++i) {
        const sym = syms[i];
        if (sym > 255) {
            const len = (sym >> 18) & 31;
            wbits16(out, p, lm[len + 257]), p += ll[len + 257];
            if (len > 7) wbits(out, p, (sym >> 23) & 31), p += fleb[len];
            const dst = sym & 31;
            wbits16(out, p, dm[dst]), p += dl[dst];
            if (dst > 3) wbits16(out, p, (sym >> 5) & 8191), p += fdeb[dst];
        } else {
            wbits16(out, p, lm[sym]), p += ll[sym];
        }
    }
    wbits16(out, p, lm[256]);
    return p + ll[256];
}

// deflate options (nice << 13) | chain
const deo = /*#__PURE__*/ new i32([65540, 131080, 131088, 131104, 262176, 1048704, 1048832, 2114560, 2117632]);

// empty
const et = /*#__PURE__*/new u8(0);

type DeflateState = {
    // head
    h?: Uint16Array;
    // prev
    p?: Uint16Array;
    // index
    i?: number;
    // end index
    z?: number;
    // wait index
    w?: number;
    // remainder byte info
    r?: number;
    // last chunk
    l: number;
};

// compresses data into a raw DEFLATE buffer
const dflt = (dat: Uint8Array, lvl: number, plvl: number, pre: number, post: number, st: DeflateState) => {
    const s = st.z || dat.length;
    const o = new u8(pre + s + 5 * (1 + Math.ceil(s / 7000)) + post);
    // writing to this writes to the output buffer
    const w = o.subarray(pre, o.length - post);
    const lst = st.l;
    let pos = (st.r || 0) & 7;
    if (lvl) {
        if (pos) w[0] = st.r >> 3;
        const opt = deo[lvl - 1];
        const n = opt >> 13, c = opt & 8191;
        const msk = (1 << plvl) - 1;
        //    prev 2-byte val map    curr 2-byte val map
        const prev = st.p || new u16(32768), head = st.h || new u16(msk + 1);
        const bs1 = Math.ceil(plvl / 3), bs2 = 2 * bs1;
        const hsh = (i: number) => (dat[i] ^ (dat[i + 1] << bs1) ^ (dat[i + 2] << bs2)) & msk;
        // 24576 is an arbitrary number of maximum symbols per block
        // 424 buffer for last block
        const syms = new i32(25000);
        // length/literal freq   distance freq
        const lf = new u16(288), df = new u16(32);
        //  l/lcnt  exbits  index          l/lind  waitdx          blkpos
        let lc = 0, eb = 0, i = st.i || 0, li = 0, wi = st.w || 0, bs = 0;
        for (; i + 2 < s; ++i) {
            // hash value
            const hv = hsh(i);
            // index mod 32768    previous index mod
            let imod = i & 32767, pimod = head[hv];
            prev[imod] = pimod;
            head[hv] = imod;
            // We always should modify head and prev, but only add symbols if
            // this data is not yet processed ("wait" for wait index)
            if (wi <= i) {
                // bytes remaining
                const rem = s - i;
                if ((lc > 7000 || li > 24576) && (rem > 423 || !lst)) {
                    pos = wblk(dat, w, 0, syms, lf, df, eb, li, bs, i - bs, pos);
                    li = lc = eb = 0, bs = i;
                    for (let j = 0; j < 286; ++j) lf[j] = 0;
                    for (let j = 0; j < 30; ++j) df[j] = 0;
                }
                //  len    dist   chain
                let l = 2, d = 0, ch = c, dif = imod - pimod & 32767;
                if (rem > 2 && hv == hsh(i - dif)) {
                    const maxn = Math.min(n, rem) - 1;
                    const maxd = Math.min(32767, i);
                    // max possible length
                    // not capped at dif because decompressors implement "rolling" index population
                    const ml = Math.min(258, rem);
                    while (dif <= maxd && --ch && imod != pimod) {
                        if (dat[i + l] == dat[i + l - dif]) {
                            let nl = 0;
                            for (; nl < ml && dat[i + nl] == dat[i + nl - dif]; ++nl);
                            if (nl > l) {
                                l = nl, d = dif;
                                // break out early when we reach "nice" (we are satisfied enough)
                                if (nl > maxn) break;
                                // now, find the rarest 2-byte sequence within this
                                // length of literals and search for that instead.
                                // Much faster than just using the start
                                const mmd = Math.min(dif, nl - 2);
                                let md = 0;
                                for (let j = 0; j < mmd; ++j) {
                                    const ti = i - dif + j & 32767;
                                    const pti = prev[ti];
                                    const cd = ti - pti & 32767;
                                    if (cd > md) md = cd, pimod = ti;
                                }
                            }
                        }
                        // check the previous match
                        imod = pimod, pimod = prev[imod];
                        dif += imod - pimod & 32767;
                    }
                }
                // d will be nonzero only when a match was found
                if (d) {
                    // store both dist and len data in one int32
                    // Make sure this is recognized as a len/dist with 28th bit (2^28)
                    syms[li++] = 268435456 | (revfl[l] << 18) | revfd[d];
                    const lin = revfl[l] & 31, din = revfd[d] & 31;
                    eb += fleb[lin] + fdeb[din];
                    ++lf[257 + lin];
                    ++df[din];
                    wi = i + l;
                    ++lc;
                } else {
                    syms[li++] = dat[i];
                    ++lf[dat[i]];
                }
            }
        }
        for (i = Math.max(i, wi); i < s; ++i) {
            syms[li++] = dat[i];
            ++lf[dat[i]];
        }
        pos = wblk(dat, w, lst, syms, lf, df, eb, li, bs, i - bs, pos);
        if (!lst) {
            st.r = (pos & 7) | w[(pos / 8) | 0] << 3;
            // shft(pos) now 1 less if pos & 7 != 0
            pos -= 7;
            st.h = head, st.p = prev, st.i = i, st.w = wi;
        }
    } else {
        for (let i = st.w || 0; i < s + lst; i += 65535) {
            // end
            let e = i + 65535;
            if (e >= s) {
                // write final block
                w[(pos / 8) | 0] = lst;
                e = s;
            }
            pos = wfblk(w, pos + 1, dat.subarray(i, e));
        }
        st.i = s;
    }
    return slc(o, 0, pre + shft(pos) + post);
}

/**
 * Options for decompressing a DEFLATE stream
 */
interface InflateStreamOptions {
    /**
     * The dictionary used to compress the original data. If no dictionary was used during compression, this option has no effect.
     * 
     * Supplying the wrong dictionary during decompression usually yields corrupt output or causes an invalid distance error.
     */
    dictionary?: Uint8Array;
}

/**
 * Options for decompressing DEFLATE data
 */
interface InflateOptions extends InflateStreamOptions {
    /**
     * The buffer into which to write the decompressed data. Saves memory if you know the decompressed size in advance.
     *
     * Note that if the decompression result is larger than the size of this buffer, it will be truncated to fit.
     */
    out?: Uint8Array;
}

/**
 * Options for compressing data into a DEFLATE format
 */
interface DeflateOptions {
    /**
     * The level of compression to use, ranging from 0-9.
     * 
     * 0 will store the data without compression.
     * 1 is fastest but compresses the worst, 9 is slowest but compresses the best.
     * The default level is 6.
     * 
     * Typically, binary data benefits much more from higher values than text data.
     * In both cases, higher values usually take disproportionately longer than the reduction in final size that results.
     * 
     * For example, a 1 MB text file could:
     * - become 1.01 MB with level 0 in 1ms
     * - become 400 kB with level 1 in 10ms
     * - become 320 kB with level 9 in 100ms
     */
    level?: 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9;
    /**
     * The memory level to use, ranging from 0-12. Increasing this increases speed and compression ratio at the cost of memory.
     * 
     * Note that this is exponential: while level 0 uses 4 kB, level 4 uses 64 kB, level 8 uses 1 MB, and level 12 uses 16 MB.
     * It is recommended not to lower the value below 4, since that tends to hurt performance.
     * In addition, values above 8 tend to help very little on most data and can even hurt performance.
     * 
     * The default value is automatically determined based on the size of the input data.
     */
    mem?: 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12;
    /**
     * A buffer containing common byte sequences in the input data that can be used to significantly improve compression ratios.
     * 
     * Dictionaries should be 32kB or smaller and include strings or byte sequences likely to appear in the input.
     * The decompressor must supply the same dictionary as the compressor to extract the original data.
     * 
     * Dictionaries only improve aggregate compression ratio when reused across multiple small inputs. They should typically not be used otherwise.
     * 
     * Avoid using dictionaries with GZIP and ZIP to maximize software compatibility.
     */
    dictionary?: Uint8Array;
};

// deflate with opts
const dopt = (dat: Uint8Array, opt: DeflateOptions, pre: number, post: number, st?: DeflateState) => {
    if (!st) {
        st = { l: 1 };
        if (opt.dictionary) {
            const dict = opt.dictionary.subarray(-32768);
            const newDat = new u8(dict.length + dat.length);
            newDat.set(dict);
            newDat.set(dat, dict.length);
            dat = newDat;
            st.w = dict.length;
        }
    }
    return dflt(dat, opt.level == null ? 6 : opt.level, opt.mem == null ? (st.l ? Math.ceil(Math.max(8, Math.min(13, Math.log(dat.length))) * 1.5) : 20) : (12 + opt.mem), pre, post, st);
}

/**
 * Compresses data with DEFLATE without any wrapper
 * @param data The data to compress
 * @param opts The compression options
 * @returns The deflated version of the data
 */
function deflateSync(data: Uint8Array, opts?: DeflateOptions) {
    return dopt(data, opts || {}, 0, 0);
}

/**
 * Expands DEFLATE data with no wrapper
 * @param data The data to decompress
 * @param opts The decompression options
 * @returns The decompressed version of the data
 */
function inflateSync(data: Uint8Array, opts?: InflateOptions) {
    return inflt(data, { i: 2 }, opts && opts.out, opts && opts.dictionary);
}