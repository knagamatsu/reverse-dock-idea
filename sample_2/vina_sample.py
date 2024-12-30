from vina import Vina

# Vinaオブジェクトの作成
v = Vina(sf_name='vina')

# レセプターとリガンドの設定
v.set_receptor('1iep_receptor.pdbqt')
v.set_ligand_from_file('1iep_ligand.pdbqt')

# グリッドボックスの設定
center = [15.190, 53.903, 16.917]
box_size = [20, 20, 20]
v.compute_vina_maps(center=center, box_size=box_size)

# ドッキングの実行
v.dock(exhaustiveness=8, n_poses=20)
v.write_poses('docked_ligands.pdbqt', n_poses=5)
