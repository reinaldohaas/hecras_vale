import re
import json
from pathlib import Path

content_path = Path(r"C:\Users\haas\.gemini\antigravity\brain\cfff4a06-4a68-4caf-8a62-fe21ae9bf18f\.system_generated\steps\2448\content.md")
with open(content_path, 'r', encoding='utf-8') as f:
    html = f.read()

pattern = re.compile(r'<tr>\s*<td[^>]*>(\d{4})</td>\s*<td[^>]*>([\d/]+)</td>\s*<td[^>]*>([\d,]+)</td>\s*</tr>')
matches = pattern.findall(html)

records = []
for ano, data, cota in matches:
    cota_float = float(cota.replace(',', '.'))
    records.append({
        'ano': int(ano),
        'data': data.strip(),
        'cota_m': cota_float
    })

records_sorted = sorted(records, key=lambda x: x['cota_m'], reverse=True)

out_file = Path(__file__).resolve().parent.parent.parent / "data" / "blumenau_103_enchentes.json"
out_file.parent.mkdir(parents=True, exist_ok=True)
with open(out_file, 'w', encoding='utf-8') as f:
    json.dump(records, f, indent=2, ensure_ascii=False)

print(f"Total de enchentes oficiais extraídas: {len(records)}")
print("\nTop 15 Maiores Enchentes da História de Blumenau (Defesa Civil):")
for i, r in enumerate(records_sorted[:15], 1):
    a = r['ano']
    d = r['data']
    c = r['cota_m']
    print(f"   {i:2d}. {a} ({d}) -> Cota: {c:5.2f} m")
