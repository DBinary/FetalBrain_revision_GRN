from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import argparse, subprocess, time, hashlib, re, json

p = argparse.ArgumentParser()
p.add_argument('url')
p.add_argument('size', type=int)
p.add_argument('publisher_bsd_sum', type=int)
p.add_argument('output')
p.add_argument('--seconds', type=int, default=900)
a = p.parse_args()
out = Path(a.output)
out.parent.mkdir(parents=True, exist_ok=True)
parts = Path('/tmp/t577_macaque_download_parts') / out.name
parts.mkdir(parents=True, exist_ok=True)
block = 1024 * 1024
total = (a.size + block - 1) // block
deadline = time.monotonic() + a.seconds
started = time.monotonic()

def get_part(i):
    start = i * block
    end = min(a.size - 1, start + block - 1)
    expected = end - start + 1
    path = parts / f'{i:06d}.bin'
    headers = parts / f'{i:06d}.headers'
    proof = parts / f'{i:06d}.subranges.json'
    content_range = f'bytes {start}-{end}/{a.size}'
    if path.is_file() and path.stat().st_size == expected and proof.is_file():
        record = json.loads(proof.read_text())
        if record['start'] == start and record['end'] == end and hashlib.sha256(path.read_bytes()).hexdigest() == record['sha256']:
            return expected
    if path.is_file() and path.stat().st_size == expected and headers.is_file():
        found = re.findall(r'(?im)^content-range:\s*([^\r\n]+)', headers.read_text())
        if found and found[-1].strip() == content_range:
            return expected
    error = ''
    for attempt in range(4):
        if time.monotonic() > deadline:
            raise TimeoutError('Range download time budget exhausted')
        partial = path.with_suffix('.partial')
        command = ['curl', '--silent', '--show-error', '--fail', '--location', '--connect-timeout', '15', '--max-time', '60', '--range', f'{start}-{end}', '--max-filesize', str(expected), '--dump-header', str(headers), '--output', str(partial), a.url]
        result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=70)
        found = re.findall(r'(?im)^content-range:\s*([^\r\n]+)', headers.read_text()) if headers.exists() else []
        valid = result.returncode == 0 and partial.is_file() and partial.stat().st_size == expected and found and found[-1].strip() == content_range
        if valid:
            partial.replace(path)
            return expected
        error = result.stderr[-300:]
        time.sleep(1)
    if expected > 65536:
        spans = []
        for sub_start in range(start, end + 1, 65536):
            sub_end = min(end, sub_start + 65535)
            sub_size = sub_end - sub_start + 1
            tiny = parts / f'{i:06d}.{sub_start}.tiny'
            tiny_headers = tiny.with_suffix('.headers')
            for attempt in range(4):
                if time.monotonic() > deadline:
                    raise TimeoutError('Subrange repair time budget exhausted')
                cmd = command.copy()
                cmd[cmd.index('--range') + 1] = f'{sub_start}-{sub_end}'
                cmd[cmd.index('--max-filesize') + 1] = str(sub_size)
                cmd[cmd.index('--dump-header') + 1] = str(tiny_headers)
                cmd[cmd.index('--output') + 1] = str(tiny)
                sub_result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=70)
                found = re.findall(r'(?im)^content-range:\s*([^\r\n]+)', tiny_headers.read_text()) if tiny_headers.exists() else []
                if sub_result.returncode == 0 and tiny.is_file() and tiny.stat().st_size == sub_size and found and found[-1].strip() == f'bytes {sub_start}-{sub_end}/{a.size}':
                    spans.append(tiny)
                    break
            else:
                raise RuntimeError(f'Failed subrange {sub_start}-{sub_end}')
        with partial.open('wb') as handle:
            for tiny in spans:
                handle.write(tiny.read_bytes())
        if partial.stat().st_size != expected:
            raise RuntimeError('Combined subrange length mismatch')
        partial.replace(path)
        proof.write_text(json.dumps({'start': start, 'end': end, 'source_url': a.url, 'validated_subranges': len(spans), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}))
        print('DOWNLOAD_PART_REPAIRED', json.dumps({'part': i, 'subranges': len(spans), 'bytes': expected}), flush=True)
        return expected
    raise RuntimeError(f'Part {i} failed after 4 attempts: {error}')

print('DOWNLOAD_START', json.dumps({'url': a.url, 'bytes': a.size, 'parts': total, 'threads': 8, 'max_seconds': a.seconds}), flush=True)
completed = 0
received = 0
failed = []
with ThreadPoolExecutor(max_workers=8) as executor:
    futures = [executor.submit(get_part, i) for i in range(total)]
    for future in as_completed(futures):
        try:
            received += future.result()
        except Exception as exc:
            failed.append(str(exc))
            print('DOWNLOAD_PART_FAILED', json.dumps({'error': str(exc)}), flush=True)
            continue
        completed += 1
        if completed % 32 == 0 or completed == total:
            print('DOWNLOAD_PROGRESS', json.dumps({'parts': completed, 'total_parts': total, 'bytes': received, 'seconds': round(time.monotonic() - started, 1)}), flush=True)
if failed:
    raise RuntimeError(json.dumps({'failed_parts': len(failed), 'details': failed[:10]}))

assembled = out.with_name(out.name + '.assembled')
hash_value = hashlib.sha256()
with assembled.open('wb') as handle:
    for i in range(total):
        chunk = (parts / f'{i:06d}.bin').read_bytes()
        hash_value.update(chunk)
        handle.write(chunk)
actual = hash_value.hexdigest()
checksum_fields = subprocess.check_output(['sum', str(assembled)], text=True).split()
bsd_sum, publisher_blocks = int(checksum_fields[0]), int(checksum_fields[1])
if assembled.stat().st_size != a.size or bsd_sum != a.publisher_bsd_sum or publisher_blocks != (a.size + 1023) // 1024:
    raise RuntimeError(f'Official Ensembl checksum or length mismatch: BSD {bsd_sum}, expected {a.publisher_bsd_sum}')
import gzip
with gzip.open(assembled, 'rb') as gzip_handle:
    for gzip_block in iter(lambda: gzip_handle.read(8*1024*1024), b''):
        pass
print('OFFICIAL_CHECKSUM_AND_GZIP_CRC_VERIFIED', json.dumps({'bsd_sum': bsd_sum, 'publisher_blocks': publisher_blocks, 'sha256_local': actual}), flush=True)
if out.exists():
    raise FileExistsError(out)
assembled.replace(out)
print('DOWNLOAD_VERIFIED', json.dumps({'output': str(out), 'bytes': out.stat().st_size, 'sha256': actual, 'seconds': round(time.monotonic() - started, 1)}), flush=True)
