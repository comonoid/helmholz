#!/usr/bin/env python3
# §888: PPM(P6) -> PNG (stdlib zlib), чтобы демонстрационные кадры открывались
# везде. Использование: ppm2png.py in.ppm out.png
import sys, zlib, struct

def read_ppm(p):
    d = open(p, 'rb').read()
    parts = d.split(b'\n', 3)
    assert parts[0] == b'P6'
    w, h = map(int, parts[1].split())
    return w, h, parts[3]

def chunk(tag, data):
    c = tag + data
    return struct.pack('>I', len(data)) + c + struct.pack('>I', zlib.crc32(c) & 0xffffffff)

def write_png(p, w, h, px):
    raw = b''.join(b'\x00' + px[y*w*3:(y+1)*w*3] for y in range(h))
    open(p, 'wb').write(
        b'\x89PNG\r\n\x1a\n' +
        chunk(b'IHDR', struct.pack('>IIBBBBB', w, h, 8, 2, 0, 0, 0)) +
        chunk(b'IDAT', zlib.compress(raw, 6)) +
        chunk(b'IEND', b''))

w, h, px = read_ppm(sys.argv[1])
write_png(sys.argv[2], w, h, px)
print(f"{sys.argv[1]} -> {sys.argv[2]} ({w}x{h})")
