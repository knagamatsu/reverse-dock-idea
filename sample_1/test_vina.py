from vina import Vina
import os
import sys

def test_vina_installation():
    try:
        # Vinaオブジェクトの作成
        v = Vina()
        version = v.version()
        print(f"Vina version: {version}")
        
        # コマンドラインツールのテスト
        stream = os.popen('vina --help')
        output = stream.read()
        if 'AutoDock Vina' in output:
            print("Command line Vina is properly installed")
        
        print("All vina tests passed successfully!")
        return True
    except Exception as e:
        print(f"Error during vina test: {str(e)}")
        return False

if __name__ == "__main__":
    success = test_vina_installation()
    sys.exit(0 if success else 1)