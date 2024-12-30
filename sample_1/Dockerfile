FROM continuumio/miniconda3

# システムの依存関係をインストール
RUN apt-get update && apt-get install -y \
    build-essential \
    libboost-all-dev \
    swig \
    wget \
    git \
    && rm -rf /var/lib/apt/lists/*

# 作業ディレクトリを設定
WORKDIR /app

# conda環境を設定
RUN conda create -n vina_env python=3.10 -y && \
    echo "conda activate vina_env" >> ~/.bashrc

# conda-forgeチャンネルを追加とパッケージのインストール
RUN . /opt/conda/etc/profile.d/conda.sh && \
    conda activate vina_env && \
    conda config --add channels conda-forge && \
    conda install -y -c conda-forge \
    rdkit \
    numpy \
    swig \
    boost-cpp \
    libboost \
    streamlit \
    pip

# Pythonパッケージのインストール
RUN . /opt/conda/etc/profile.d/conda.sh && \
    conda activate vina_env && \
    pip install --no-cache-dir \
    streamlit-ketcher \
    py3dmol \
    meeko \
    vina

# AutoDock Vinaバイナリをダウンロードしてインストール
RUN wget https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.3/vina_1.2.3_linux_x86_64 && \
    chmod +x vina_1.2.3_linux_x86_64 && \
    mv vina_1.2.3_linux_x86_64 /usr/local/bin/vina

# ポートを公開
EXPOSE 8501

# アプリケーションファイルをコピー
COPY . .

# Streamlitを起動
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "vina_env", "streamlit", "run", "app.py", "--server.address", "0.0.0.0"]