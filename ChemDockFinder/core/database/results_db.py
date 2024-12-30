from typing import Dict, List, Optional
import pathlib
import pandas as pd
import json
import logging
from datetime import datetime

logger = logging.getLogger(__name__)

class ResultsDatabase:
    """ドッキング結果のデータベース管理クラス"""
    
    def __init__(self, base_dir: Optional[pathlib.Path] = None):
        self.base_dir = base_dir or pathlib.Path('./data/results')
        self.results_dir = self.base_dir / 'docking'
        self.metadata_file = self.base_dir / 'results_metadata.json'
        
        # ディレクトリの初期化
        self.results_dir.mkdir(parents=True, exist_ok=True)
        
        # メタデータの読み込み
        self.metadata = self._load_metadata()

    def _load_metadata(self) -> Dict:
        """メタデータの読み込み"""
        if self.metadata_file.exists():
            with open(self.metadata_file, 'r') as f:
                return json.load(f)
        return {'calculations': {}}

    def _save_metadata(self):
        """メタデータの保存"""
        with open(self.metadata_file, 'w') as f:
            json.dump(self.metadata, f, indent=2)

    def add_calculation_result(self, 
                             calculation_id: str,
                             results: Dict) -> bool:
        """
        計算結果の追加
        
        Args:
            calculation_id: 計算ID
            results: 計算結果
        """
        try:
            # 結果ディレクトリの作成
            result_dir = self.results_dir / calculation_id
            result_dir.mkdir(exist_ok=True)
            
            # 結果ファイルの保存
            result_file = result_dir / 'results.json'
            with open(result_file, 'w') as f:
                json.dump(results, f, indent=2)
            
            # メタデータの更新
            self.metadata['calculations'][calculation_id] = {
                'date': datetime.now().isoformat(),
                'smiles': results.get('smiles', ''),
                'status': results.get('status', 'completed'),
                'result_file': str(result_file)
            }
            
            self._save_metadata()
            return True
            
        except Exception as e:
            logger.error(f"Error adding calculation result {calculation_id}: {e}")
            return False

    def get_calculation_result(self, calculation_id: str) -> Optional[Dict]:
        """計算結果の取得"""
        try:
            meta = self.metadata['calculations'].get(calculation_id)
            if meta and pathlib.Path(meta['result_file']).exists():
                with open(meta['result_file'], 'r') as f:
                    return json.load(f)
            return None
        except Exception as e:
            logger.error(f"Error getting calculation result {calculation_id}: {e}")
            return None

    def get_all_calculations(self) -> pd.DataFrame:
        """全計算結果の概要を取得"""
        data = []
        for calc_id, info in self.metadata['calculations'].items():
            row = {
                'Calculation ID': calc_id,
                'Date': info['date'],
                'SMILES': info['smiles'],
                'Status': info['status']
            }
            data.append(row)
        return pd.DataFrame(data)

    def get_recent_calculations(self, limit: int = 10) -> pd.DataFrame:
        """最近の計算結果を取得"""
        df = self.get_all_calculations()
        if not df.empty:
            df['Date'] = pd.to_datetime(df['Date'])
            return df.sort_values('Date', ascending=False).head(limit)
        return df

    def remove_calculation(self, calculation_id: str) -> bool:
        """計算結果の削除"""
        try:
            if calculation_id in self.metadata['calculations']:
                result_dir = self.results_dir / calculation_id
                if result_dir.exists():
                    import shutil
                    shutil.rmtree(result_dir)
                del self.metadata['calculations'][calculation_id]
                self._save_metadata()
                return True
            return False
        except Exception as e:
            logger.error(f"Error removing calculation {calculation_id}: {e}")
            return False

    def get_calculation_summary(self, calculation_id: str) -> Optional[Dict]:
        """計算結果の要約を取得"""
        try:
            result = self.get_calculation_result(calculation_id)
            if result:
                if 'results' in result:
                    energies = [r.get('best_energy', 0) for r in result['results'] 
                              if 'error' not in r]
                    return {
                        'id': calculation_id,
                        'date': self.metadata['calculations'][calculation_id]['date'],
                        'smiles': result.get('smiles', ''),
                        'num_proteins': len(result.get('results', [])),
                        'best_energy': min(energies) if energies else None,
                        'status': result.get('status', 'unknown')
                    }
            return None
        except Exception as e:
            logger.error(f"Error getting calculation summary {calculation_id}: {e}")
            return None