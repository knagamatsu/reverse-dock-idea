from typing import Dict, List, Optional, Union
import pathlib
import pandas as pd
import json
import logging
from datetime import datetime
from .vina import VinaManager

# ロガーの設定
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class DockingCalculator:
    """ドッキング計算の制御クラス"""
    
    def __init__(self, base_dir: Optional[pathlib.Path] = None):
        """
        初期化
        
        Args:
            base_dir: 結果保存用のベースディレクトリ
        """
        self.vina = VinaManager()
        self.base_dir = base_dir or pathlib.Path('./data/results')
        self.results_dir = self.base_dir / 'docking'
        self.logs_dir = self.base_dir / 'logs'
        
        # ディレクトリの作成
        self.results_dir.mkdir(parents=True, exist_ok=True)
        self.logs_dir.mkdir(parents=True, exist_ok=True)

    def _save_calculation_log(self, 
                            calculation_id: str, 
                            data: Dict) -> pathlib.Path:
        """計算ログの保存"""
        log_file = self.logs_dir / f"{calculation_id}.json"
        with open(log_file, 'w') as f:
            json.dump(data, f, indent=2)
        return log_file

    def prepare_calculation(self, 
                          smiles: str, 
                          protein_files: Union[pathlib.Path, List[pathlib.Path]], 
                          params: Dict) -> str:
        """
        ドッキング計算の準備
        
        Args:
            smiles: リガンドのSMILES
            protein_files: タンパク質構造ファイル（単一または複数）
            params: 計算パラメータ
            
        Returns:
            str: 計算ID
        """
        try:
            # 単一ファイルの場合はリストに変換
            if isinstance(protein_files, pathlib.Path):
                protein_files = [protein_files]

            # 計算IDの生成
            calculation_id = datetime.now().strftime('%Y%m%d_%H%M%S')
            
            # 計算情報の保存
            calculation_info = {
                'id': calculation_id,
                'smiles': smiles,
                'proteins': [str(p) for p in protein_files],
                'params': params,
                'status': 'prepared',
                'timestamp': datetime.now().isoformat()
            }
            
            self._save_calculation_log(calculation_id, calculation_info)
            
            return calculation_id
            
        except Exception as e:
            logger.error(f"Error in prepare_calculation: {e}")
            raise

    def run_docking(self, 
                   calculation_id: str,
                   update_callback: Optional[callable] = None) -> Dict:
        """
        ドッキング計算の実行
        
        Args:
            calculation_id: 計算ID
            update_callback: 進捗更新用コールバック関数
            
        Returns:
            Dict: 計算結果
        """
        try:
            # 計算情報の読み込み
            log_file = self.logs_dir / f"{calculation_id}.json"
            with open(log_file, 'r') as f:
                calculation_info = json.load(f)
            
            if calculation_info['status'] != 'prepared':
                raise ValueError(f"Invalid calculation status: {calculation_info['status']}")
            
            # 結果保存用ディレクトリ
            calc_dir = self.results_dir / calculation_id
            calc_dir.mkdir(exist_ok=True)
            
            results = []
            total_proteins = len(calculation_info['proteins'])
            
            for i, protein_file in enumerate(calculation_info['proteins'], 1):
                try:
                    # 進捗更新
                    if update_callback:
                        progress = (i / total_proteins) * 100
                        update_callback(progress)
                    
                    # ドッキング実行
                    result_file, result_data = self.vina.run_docking(
                        protein_file=pathlib.Path(protein_file),
                        ligand_smiles=calculation_info['smiles'],
                        params=calculation_info['params']
                    )
                    
                    # 結果ファイルの移動
                    new_result_file = calc_dir / f"{pathlib.Path(protein_file).stem}_result.pdbqt"
                    result_file.rename(new_result_file)
                    
                    # 結果の集計
                    result_info = {
                        'protein': protein_file,
                        'result_file': str(new_result_file),
                        'best_energy': min(result_data['energies']),
                        'all_energies': result_data['energies'],
                        'center': result_data['center'],
                        'box_size': result_data['box_size']
                    }
                    
                    results.append(result_info)
                    
                except Exception as e:
                    logger.error(f"Error in docking with {protein_file}: {e}")
                    result_info = {
                        'protein': protein_file,
                        'error': str(e)
                    }
                    results.append(result_info)
            
            # 結果の保存
            calculation_info['status'] = 'completed'
            calculation_info['results'] = results
            calculation_info['completion_time'] = datetime.now().isoformat()
            
            self._save_calculation_log(calculation_id, calculation_info)
            
            # 結果のDataFrame作成
            df = pd.DataFrame([
                r for r in results if 'error' not in r
            ])
            
            if not df.empty:
                df.sort_values('best_energy', ascending=True, inplace=True)
                csv_file = calc_dir / 'results.csv'
                df.to_csv(csv_file, index=False)
            
            return calculation_info
            
        except Exception as e:
            logger.error(f"Error in run_docking: {e}")
            raise

    def get_calculation_status(self, calculation_id: str) -> Dict:
        """計算状態の取得"""
        try:
            log_file = self.logs_dir / f"{calculation_id}.json"
            with open(log_file, 'r') as f:
                return json.load(f)
        except Exception as e:
            logger.error(f"Error in get_calculation_status: {e}")
            raise

    def get_results_dataframe(self, calculation_id: str) -> Optional[pd.DataFrame]:
        """結果のDataFrame取得"""
        try:
            csv_file = self.results_dir / calculation_id / 'results.csv'
            if csv_file.exists():
                return pd.read_csv(csv_file)
            return None
        except Exception as e:
            logger.error(f"Error in get_results_dataframe: {e}")
            raise

    def cleanup(self):
        """クリーンアップ処理"""
        try:
            self.vina.cleanup()
        except Exception as e:
            logger.error(f"Error in cleanup: {e}")
            raise