from typing import Dict, List, Tuple, Optional
import numpy as np
from vina import Vina
from rdkit import Chem
from rdkit.Chem import AllChem
import tempfile
import pathlib
import logging
from datetime import datetime

# ロガーの設定
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class VinaManager:
    """AutoDock Vinaの管理クラス"""
    
    def __init__(self):
        """初期化"""
        self.vina = Vina(sf_name='vina')
        self.temp_dir = pathlib.Path(tempfile.mkdtemp())
        logger.info(f"Created temporary directory: {self.temp_dir}")

    def prepare_ligand(self, smiles: str) -> pathlib.Path:
        """
        SMILESからリガンドのPDBQT作成
        
        Args:
            smiles (str): 入力SMILES
            
        Returns:
            pathlib.Path: 生成されたPDBQTファイルのパス
            
        Raises:
            Exception: 変換処理でエラーが発生した場合
        """
        try:
            # RDKitで3D構造生成
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError(f"Invalid SMILES: {smiles}")
                
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, randomSeed=42)
            AllChem.MMFFOptimizeMolecule(mol)
            
            # 一時ファイルの作成
            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
            ligand_pdb = self.temp_dir / f"ligand_{timestamp}.pdb"
            ligand_pdbqt = self.temp_dir / f"ligand_{timestamp}.pdbqt"
            
            # PDBファイルの作成
            with open(ligand_pdb, 'w') as f:
                f.write(Chem.MolToPDBBlock(mol))
            
            # PDBQTへの変換（meekoを使用）
            self.vina.convert_pdb_to_pdbqt(str(ligand_pdb), str(ligand_pdbqt))
            
            return ligand_pdbqt
        
        except Exception as e:
            logger.error(f"Error in prepare_ligand: {e}")
            raise

    def prepare_protein(self, pdb_file: pathlib.Path) -> pathlib.Path:
        """
        タンパク質構造のPDBQT作成
        
        Args:
            pdb_file (pathlib.Path): 入力PDBファイルのパス
            
        Returns:
            pathlib.Path: 生成されたPDBQTファイルのパス
            
        Raises:
            Exception: 変換処理でエラーが発生した場合
        """
        try:
            if not pdb_file.exists():
                raise FileNotFoundError(f"Protein PDB file not found: {pdb_file}")
                
            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
            protein_pdbqt = self.temp_dir / f"protein_{timestamp}.pdbqt"
            
            # PDBQTへの変換
            self.vina.convert_pdb_to_pdbqt(str(pdb_file), str(protein_pdbqt))
            
            return protein_pdbqt
        
        except Exception as e:
            logger.error(f"Error in prepare_protein: {e}")
            raise

    def calculate_box_center_size(self, protein_file: pathlib.Path) -> Tuple[List[float], List[float]]:
        """
        タンパク質構造から計算ボックスの中心と大きさを計算
        
        Args:
            protein_file (pathlib.Path): タンパク質構造ファイルのパス
            
        Returns:
            Tuple[List[float], List[float]]: (中心座標, ボックスサイズ)
            
        Raises:
            Exception: 計算でエラーが発生した場合
        """
        try:
            coords = []
            with open(protein_file, 'r') as f:
                for line in f:
                    if line.startswith(('ATOM', 'HETATM')):
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        coords.append([x, y, z])
            
            if not coords:
                raise ValueError(f"No atom coordinates found in {protein_file}")
            
            coords = np.array(coords)
            center = coords.mean(axis=0).tolist()
            
            # タンパク質全体を覆うサイズ + パディング(10Å)
            size = ((coords.max(axis=0) - coords.min(axis=0)) + 10.0).tolist()
            
            return center, size
        
        except Exception as e:
            logger.error(f"Error in calculate_box_center_size: {e}")
            raise

    def run_docking(self, 
                   protein_file: pathlib.Path,
                   ligand_smiles: str,
                   params: Dict) -> Tuple[pathlib.Path, Dict]:
        """
        ドッキング計算の実行
        
        Args:
            protein_file (pathlib.Path): タンパク質構造ファイルのパス
            ligand_smiles (str): リガンドのSMILES
            params (Dict): 計算パラメータ
            
        Returns:
            Tuple[pathlib.Path, Dict]: (結果ファイルのパス, 計算結果の辞書)
            
        Raises:
            Exception: 計算でエラーが発生した場合
        """
        try:
            # リガンドとタンパク質の準備
            ligand_pdbqt = self.prepare_ligand(ligand_smiles)
            protein_pdbqt = self.prepare_protein(protein_file)
            
            # 計算ボックスの設定
            center, box_size = self.calculate_box_center_size(protein_file)
            
            # パラメータの取得
            exhaustiveness = params.get('exhaustiveness', 8)
            num_modes = params.get('num_modes', 9)
            energy_range = params.get('energy_range', 3)
            
            # 指定されたグリッドサイズがある場合は使用
            if 'grid_size' in params:
                box_size = [params['grid_size']] * 3
            
            # Vinaの設定
            self.vina.set_receptor(str(protein_pdbqt))
            self.vina.set_ligand_from_file(str(ligand_pdbqt))
            
            # グリッドの計算
            self.vina.compute_vina_maps(
                center=center,
                box_size=box_size
            )
            
            # ドッキング実行
            self.vina.dock(
                exhaustiveness=exhaustiveness,
                n_poses=num_modes,
                energy_range=energy_range
            )
            
            # 結果の保存
            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
            result_file = self.temp_dir / f"result_{timestamp}.pdbqt"
            self.vina.write_poses(str(result_file))
            
            # 結果の集計
            results = {
                'center': center,
                'box_size': box_size,
                'poses': self.vina.poses(),
                'energies': [pose.energy for pose in self.vina.poses()],
                'params': {
                    'exhaustiveness': exhaustiveness,
                    'num_modes': num_modes,
                    'energy_range': energy_range
                }
            }
            
            return result_file, results
        
        except Exception as e:
            logger.error(f"Error in run_docking: {e}")
            raise

    def cleanup(self):
        """一時ファイルの削除"""
        try:
            import shutil
            shutil.rmtree(self.temp_dir)
            logger.info(f"Cleaned up temporary directory: {self.temp_dir}")
        except Exception as e:
            logger.error(f"Error in cleanup: {e}")
            raise