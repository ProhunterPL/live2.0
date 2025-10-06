LIVE 2.0 → PubChem matcher (spec & mini-README)

Poniżej gotowy plik MD do repo. Po zrzuceniu obrazka klastra (PNG/JPG) skrypt:

weźmie JSON metadane o klastrze (zapisane obok obrazka),

zbuduje SMILES,

znajdzie najbardziej podobną molekułę w PubChem,

wygeneruje obraz porównawczy: po lewej klaster z Live, po prawej PubChem.

Uwaga praktyczna: z samego obrazka ciężko odtworzyć graf. Dlatego przy zapisie obrazka zapisz równolegle *.json z topologią (węzły+wiązania). Triggerem jest zapis obrazka – skrypt oczekuje też pary *.json o tej samej nazwie.

🧭 Flow

Frontend Live 2.0

Po kliknięciu „Save cluster” zapisz:

clustershots/<timestamp>.png (wizualizacja)

clustershots/<timestamp>.json (graf; patrz „Kontrakt danych”)

Watcher (backend)

Monitoruje folder clustershots/ (nowe obrazy → start joba).

Job:

wczytuje *.json, tworzy RDKit Mol → SMILES,

PubChem similarity search → wybiera top wynik,

renderuje obraz PubChem (RDKit/PubChem PNG),

składa panel porównawczy (PIL) → matches/<timestamp>_match.png,

zwraca JSON z wynikami (matches/<timestamp>_match.json).

Frontend

Po zakończeniu joba wyświetla gotowy panel i metadane.

📦 Struktura
live-2.0-chem-matcher/
├─ clustershots/              # wejście (PNG/JPG + JSON)
├─ matches/                   # wyjście (panel + result JSON)
├─ matcher/
│  ├─ matcher.py              # główny CLI
│  ├─ watcher.py              # watcher folderu
│  ├─ chem.py                 # rdkit + pubchem utils
│  └─ compose.py              # składanie panelu
├─ requirements.txt
└─ README_MATCHER.md          # ten plik

📑 Kontrakt danych (JSON obok obrazka)

Minimalny, prosty i elastyczny:

{
  "id": "cluster_2025-10-08T21-53-02",
  "nodes": [
    {"id": 0, "label": "A", "mass": 1.0, "pos": [x, y], "energy": 59.4},
    {"id": 1, "label": "A", "mass": 1.0, "pos": [x, y], "energy": 56.2}
  ],
  "bonds": [
    {"a": 0, "b": 2, "order": 1},
    {"a": 2, "b": 5, "order": 1}
  ],
  "metadata": {
    "size": 8,
    "bonds": 10,
    "density": 0.357,
    "avg_mass": 1.0,
    "total_energy": 501.78,
    "avg_energy": 62.72
  }
}


label – pseudo-pierwiastek (np. "A", "B") albo od razu "C", "H", "O".

order – 1/2/3 (jeśli nie masz, przyjmij 1).

pos tylko informacyjne (do podpisów), nieobowiązkowe.

🧪 Instalacja

requirements.txt:

rdkit-pypi==2022.9.5
requests>=2.32
pillow>=10.4
watchdog>=4.0


RDKit z pip (rdkit-pypi) jest ok do wykorzystań bez chemoinformatycznych „ciężarów”. Jeśli wolisz conda – użyj conda install -c conda-forge rdkit.

🧰 Implementacja (skróty najważniejszych fragmentów)
matcher/chem.py
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import requests, urllib.parse

# Heurystyka mapowania gdy brak realnych symboli:
# - stopień węzła: 1→H, 2→O/N, >=3→C (fallback na C)
def choose_symbol(deg: int, given_label: str | None = None) -> str:
    if given_label and given_label in {"H","C","N","O","S","P","F","Cl","Br","I"}:
        return given_label
    if deg <= 1: return "H"
    if deg == 2: return "O"
    return "C"

def json_to_mol(cluster_json: dict) -> Chem.Mol:
    nodes = cluster_json["nodes"]
    bonds = cluster_json["bonds"]
    deg = {n["id"]: 0 for n in nodes}
    for b in bonds:
        deg[b["a"]] += 1; deg[b["b"]] += 1

    rw = Chem.RWMol()
    idx_map = {}

    # add atoms
    for n in nodes:
        sym = choose_symbol(deg[n["id"]], n.get("label"))
        idx_map[n["id"]] = rw.AddAtom(Chem.Atom(sym))

    # add bonds
    order_map = {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE}
    for b in bonds:
        rw.AddBond(idx_map[b["a"]], idx_map[b["b"]], order_map.get(b.get("order",1), Chem.BondType.SINGLE))

    mol = rw.GetMol()
    AllChem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_KEKULIZE)
    return mol

def mol_to_smiles(mol: Chem.Mol) -> str:
    return Chem.MolToSmiles(mol, canonical=True)

def pubchem_similar_top(smiles: str, threshold: int = 90):
    q = urllib.parse.quote(smiles)
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/similarity/smiles/{q}/JSON?Threshold={threshold}"
    r = requests.get(url, timeout=25)
    r.raise_for_status()
    cids = r.json().get("IdentifierList", {}).get("CID", [])
    if not cids: return None
    cid = cids[0]

    props_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/IUPACName,MolecularFormula,MolecularWeight,CanonicalSMILES,InChIKey/JSON"
    p = requests.get(props_url, timeout=25).json()["PropertyTable"]["Properties"][0]
    return {
        "cid": cid,
        "name": p.get("IUPACName"),
        "formula": p.get("MolecularFormula"),
        "mw": p.get("MolecularWeight"),
        "smiles": p.get("CanonicalSMILES"),
        "inchikey": p.get("InChIKey")
    }

def render_pubchem_png(smiles: str, out_png: str, size: int = 512):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("PubChem SMILES could not be parsed")
    Draw.MolToFile(mol, out_png, size=(size, size))

matcher/compose.py
from PIL import Image, ImageDraw, ImageFont

def compose_panel(left_img_path: str, right_img_path: str, title_left: str, title_right: str, out_path: str):
    left = Image.open(left_img_path).convert("RGBA")
    right = Image.open(right_img_path).convert("RGBA")

    H = max(left.height, right.height)
    scale_h = 640  # unifikacja
    def scale(im):
        r = scale_h / im.height
        return im.resize((int(im.width*r), scale_h), Image.LANCZOS)

    left = scale(left); right = scale(right)

    gap = 40
    W = left.width + right.width + gap*3
    canvas = Image.new("RGBA", (W, scale_h + 140), (22,22,22,255))
    draw = ImageDraw.Draw(canvas)

    # tytuły
    font = ImageFont.load_default()
    draw.text((gap, 16), title_left, (255,255,255), font=font)
    draw.text((gap + left.width + gap, 16), title_right, (255,255,255), font=font)

    # obrazy
    canvas.paste(left, (gap, 80), left)
    canvas.paste(right, (gap + left.width + gap, 80), right)

    canvas.save(out_path)

matcher/matcher.py (CLI dla jednego pliku)
import json, sys, pathlib, time
from chem import json_to_mol, mol_to_smiles, pubchem_similar_top, render_pubchem_png
from compose import compose_panel
from rdkit.Chem import Draw

def run_for_pair(img_path: str):
    p = pathlib.Path(img_path)
    meta = p.with_suffix(".json")
    assert meta.exists(), f"Missing meta JSON next to image: {meta}"

    data = json.loads(meta.read_text(encoding="utf-8"))
    mol = json_to_mol(data)
    smiles = mol_to_smiles(mol)

    # render lokalnej "molekuły" z Live (z grafu)
    tmp_left = p.parent / f"{p.stem}_live_rdkit.png"
    Draw.MolToFile(mol, str(tmp_left), size=(512,512))

    # PubChem query
    top = pubchem_similar_top(smiles, threshold=80) or {}
    pubchem_png = p.parent / f"{p.stem}_pubchem.png"
    if top.get("smiles"):
        render_pubchem_png(top["smiles"], str(pubchem_png), size=512)
    else:
        # fallback: pusta karta
        pubchem_png = tmp_left

    out_dir = pathlib.Path("matches"); out_dir.mkdir(exist_ok=True)
    panel = out_dir / f"{p.stem}_match.png"
    title_left = "LIVE 2.0 cluster"
    title_right = f"PubChem CID {top.get('cid','–')} | {top.get('name','unknown')}"
    compose_panel(str(tmp_left), str(pubchem_png), title_left, title_right, str(panel))

    # zapisz wynik JSON
    result = {
        "cluster_id": data.get("id"),
        "smiles_live": smiles,
        "pubchem_top": top,
        "panel_path": str(panel)
    }
    (out_dir / f"{p.stem}_match.json").write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False))

if __name__ == "__main__":
    run_for_pair(sys.argv[1])

matcher/watcher.py (monitorowanie folderu)
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler
import time, pathlib, subprocess, sys

WATCH_DIR = "clustershots"

class Handler(FileSystemEventHandler):
    def on_created(self, event):
        p = pathlib.Path(event.src_path)
        if p.suffix.lower() in {".png",".jpg",".jpeg"}:
            json_path = p.with_suffix(".json")
            # Poczekaj chwilę aż JSON zostanie zapisany
            for _ in range(40):
                if json_path.exists(): break
                time.sleep(0.1)
            try:
                subprocess.run([sys.executable, "matcher/matcher.py", str(p)], check=True)
            except Exception as e:
                print("MATCH ERROR:", e)

if __name__ == "__main__":
    pathlib.Path(WATCH_DIR).mkdir(exist_ok=True)
    ob = Observer()
    ob.schedule(Handler(), WATCH_DIR, recursive=False)
    ob.start()
    print(f"Watching {WATCH_DIR} ...")
    try:
        while True: time.sleep(1)
    except KeyboardInterrupt:
        ob.stop(); ob.join()

⚙️ Użycie

Zainstaluj zależności

pip install -r requirements.txt


Uruchom watcher

python matcher/watcher.py


Zapisz z frontendu (przykład):

obraz: clustershots/2025-10-08_215302.png

json: clustershots/2025-10-08_215302.json

Efekt

panel: matches/2025-10-08_215302_match.png

dane: matches/2025-10-08_215302_match.json

🖥️ Hook po stronie frontendu (pseudo-TS)

Po udanym zapisie canvasa:

await savePng(`clustershots/${stamp}.png`)
await saveJson(`clustershots/${stamp}.json`, clusterGraph) // nodes+bonds
// watcher odpali matcher i zapisze wynik do /matches


Możesz też dodać polling w UI, który pokaże matches/${stamp}_match.png gdy się pojawi.

🧠 Heurystyka mapowania → „chemia”

Jeśli Twoje węzły nie mają realnych symboli pierwiastków:

stopień 1 → „H”;

stopień 2 → „O” (lub „N” – możesz dodać losowanie/proporcje);

stopień ≥3 → „C”.
To tylko proxy, ale wystarczy, żeby similarity w PubChem zwracało sensowne, proste szkielety.

Chcesz lepiej? Dodaj:

reguły wg energii/massy (np. wysokie energie → heteroatomy: O/N/S),

ręczne mapy (label "6.8" → "C", "5.2" → "O", itp.).

🔗 PubChem PUG REST (używane endpointy)

Similarity by SMILES → CID:
.../rest/pug/compound/similarity/smiles/<SMILES>/JSON?Threshold=80

Properties by CID:
.../rest/pug/compound/cid/<CID>/property/IUPACName,MolecularFormula,MolecularWeight,CanonicalSMILES,InChIKey/JSON

✅ Checklist wdrożenia

 Front zapisuje PNG + JSON (ta sama nazwa).

 Watcher działa w tle (systemd/PM2/Docker).

 RDKit zainstalowany, requirements.txt przechodzi.

 Foldery clustershots/ i matches/ istnieją (uprawnienia zapisu).

 UI nasłuchuje pojawienia się *_match.png i pokazuje użytkownikowi.

📎 Notatki końcowe

Jeśli chcesz twardą identyczność zamiast podobieństwa: po wygenerowaniu SMILES zrób zapytanie po InChIKey albo exact smiles match w PubChem.

W przyszłości możesz oceniać Tanimoto lokalnie (RDKit) na dużych setach referencyjnych i dopiero top-k sprawdzać w PubChem.

Dla materiałów okresowych (kratki) użyj innego pipeline’u (pymatgen + Materials Project) – tutaj celowo skupiamy się na klastrach „molekularnych”.

Gotowe. Wrzucasz to do repo i masz automatyczny „what does this look like in real chemistry?” dla swoich klastrów.