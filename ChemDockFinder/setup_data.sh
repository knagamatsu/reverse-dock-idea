#!/bin/bash
# setup_data.sh

# ディレクトリ作成
mkdir -p data/proteins/original
mkdir -p data/proteins/prepared
mkdir -p data/results/docking
mkdir -p data/results/logs
mkdir -p data/configs

# サンプルPDBファイルのダウンロード
cd data/proteins/original

# タンパク質構造のダウンロード
wget https://files.rcsb.org/download/1iep.pdb
wget https://files.rcsb.org/download/4gid.pdb
wget https://files.rcsb.org/download/6lu7.pdb
wget https://files.rcsb.org/download/6m0j.pdb
wget https://files.rcsb.org/download/7bz5.pdb

# サンプルリガンドのダウンロード
wget https://raw.githubusercontent.com/ccsb-scripps/AutoDock-Vina/develop/example/python_scripting/1iep_ligand.pdbqt

# メタデータファイルの作成
cd ../../..
cat > data/proteins/metadata.json << EOL
{
  "proteins": {
    "1iep": {
      "file": "data/proteins/original/1iep.pdb",
      "info": {
        "Resolution": 1.8,
        "Method": "X-ray",
        "Species": "H. sapiens",
        "Ligand": "NAD"
      }
    },
    "4gid": {
      "file": "data/proteins/original/4gid.pdb",
      "info": {
        "Resolution": 2.0,
        "Method": "X-ray",
        "Species": "H. sapiens",
        "Ligand": "ATP"
      }
    },
    "6lu7": {
      "file": "data/proteins/original/6lu7.pdb",
      "info": {
        "Resolution": 2.1,
        "Method": "X-ray",
        "Species": "SARS-CoV-2",
        "Ligand": "N3"
      }
    },
    "6m0j": {
      "file": "data/proteins/original/6m0j.pdb",
      "info": {
        "Resolution": 2.45,
        "Method": "X-ray",
        "Species": "SARS-CoV-2",
        "Ligand": "ACE2"
      }
    },
    "7bz5": {
      "file": "data/proteins/original/7bz5.pdb",
      "info": {
        "Resolution": 1.95,
        "Method": "X-ray",
        "Species": "H. sapiens",
        "Ligand": "REM"
      }
    }
  }
}
EOL

# サンプル設定ファイルの作成
cat > data/configs/default_settings.yaml << EOL
docking:
  grid_size: 24
  exhaustiveness: 8
  num_modes: 9
  energy_range: 3
  center:
    x: 0.0
    y: 0.0
    z: 0.0
  scoring: "Vina"
  use_gpu: false
  consider_waters: false
  cpu_threads: 4
EOL

echo "Data setup completed!"