from typing import Dict, List, Optional
import pathlib
import pandas as pd
import json
import logging
from datetime import datetime

logger = logging.getLogger(__name__)

class ProteinDatabase:
    """タンパク質構造データベースの管理クラス"""
    
    def __init__(self, base_dir: Optional[pathlib.Path] = None):
        self.base_dir = base_dir or pathlib.Path('./data/proteins')
        self.original_dir = self.base_dir / 'original'
        self.prepared_dir = self.base_dir / 'prepared'
        self.metadata_file = self.base_dir / 'metadata.json'
        
        # ディレクトリの初期化
        self.original_dir.mkdir(parents=True, exist_ok=True)
        self.prepared_dir.mkdir(parents=True, exist_ok=True)
        
        # メタデータの読み込み
        self.metadata = self._load_metadata()

    def _load_metadata(self) -> Dict:
        """メタデータの読み込み"""
        if self.metadata_file.exists():
            with open(self.metadata_file, 'r') as f:
                return json.load(f)
        return {'proteins': {}}

    def _save_metadata(self):
        """メタデータの保存"""
        with open(self.metadata_file, 'w') as f:
            json.dump(self.metadata, f, indent=2)

    def add_protein(self, 
                   pdb_id: str, 
                   pdb_file: pathlib.Path,
                   info: Dict) -> bool:
        """
        タンパク質構造の追加
        
        Args:
            pdb_id: PDB ID
            pdb_file: PDBファイルパス
            info: タンパク質の情報
        """
        try:
            # ファイルのコピー
            dest_file = self.original_dir / f"{pdb_id}.pdb"
            with open(pdb_file, 'rb') as src, open(dest_file, 'wb') as dst:
                dst.write(src.read())
            
            # メタデータの更新
            self.metadata['proteins'][pdb_id] = {
                'file': str(dest_file),
                'added_date': datetime.now().isoformat(),
                'info': info
            }
            
            self._save_metadata()
            return True
            
        except Exception as e:
            logger.error(f"Error adding protein {pdb_id}: {e}")
            return False

    def get_protein_info(self, pdb_id: str) -> Optional[Dict]:
        """タンパク質情報の取得"""
        return self.metadata['proteins'].get(pdb_id)

    def get_all_proteins(self) -> pd.DataFrame:
        """全タンパク質情報の取得"""
        data = []
        for pdb_id, info in self.metadata['proteins'].items():
            row = {'PDB ID': pdb_id}
            row.update(info['info'])
            data.append(row)
        return pd.DataFrame(data)

    def get_protein_file(self, pdb_id: str) -> Optional[pathlib.Path]:
        """タンパク質構造ファイルのパスを取得"""
        info = self.get_protein_info(pdb_id)
        if info:
            return pathlib.Path(info['file'])
        return None

    def remove_protein(self, pdb_id: str) -> bool:
        """タンパク質構造の削除"""
        try:
            if pdb_id in self.metadata['proteins']:
                file_path = pathlib.Path(self.metadata['proteins'][pdb_id]['file'])
                if file_path.exists():
                    file_path.unlink()
                del self.metadata['proteins'][pdb_id]
                self._save_metadata()
                return True
            return False
        except Exception as e:
            logger.error(f"Error removing protein {pdb_id}: {e}")
            return False