# Passo 1: Iniciar o Conda
eval "$(/home/lucas/anaconda3/bin/conda shell.bash hook)"

# Passo 2: Criar diretório chamado "Mapeamento".
mkdir Mapeamento

# Passo 3: Fixar o diretório "Mapeamento" como diretório de tranalho.
cd Mapeamento

# Passo 4: Criar diretório chamado "Dados" para armazenar os dados.
mkdir Dados

# Passo 5: Criar subdiretório "Arquivos Para Analise" para armazenar os dados em fastq. [SALVAR AQUI DADOS A SEREM ANALISADOS EM FORMATO FASTQ- .fastqsanger]
mkdir Dados/ArquivosParaAnalise

# Passo 6: Criar subdiretório "Sequencia de Referencia" para armazenar a sequencia de referência, fasta e gtf. [SALVAR AQUI A SEQUENCIA DE REFERÊNCIA EM FORMATO FASTA-.fa E CORRESPONDENTE GTF -.fai]
mkdir Dados/SequenciadeReferencia

# Passo 7: Verificar qualidade das sequências com FastQC. O FastQC atuará em todos os documentos no diretório. Se der certo, haverá um arquivo html para cada arquivo em fastq no diretório. Se o FastQC não tiver sido instalado, isso será informado pelo terminal, bem como o comando para instalá-lo.
fastqc Dados/ArquivosParaAnalise/*

# Passo 8: Agregar resultados de todos os documentos no diretório em apenas um documento html,usando o Multiqc. Esse passo visa a facilitar a leitura de muitos arquivos em simultâneo. Se der certo, surgirá um arquivo em html no diretório>
multiqc Dados/ArquivosParaAnalise/*fastqc*

# Passo 9: Criar diretório "Dados do Star",para armazenar os arquivos de mapeamento do STAR.
mkdir Dados/DadosdoStar

# Passo 10: Abrir o GNU para criar ídice de genomas para poder mapear reads com STAR.
nano generate_DadosdoStar.sh

# Passo 11: Copiar o código a seguir no GNU para poder criar ídice de genomas para poder mapear reads com STAR. Para salvar, clicar ctrl+o depois enter e ctrl+x para fechar o GNU.
#!/bin/bash

STAR --runThreadN 4 \
     --limitGenomeGenerateRAM 2000000000 \
     --runMode genomeGenerate \
     --genomeDir Dados/DadosdoStar \
     --genomeFastaFiles Dados/SequenciadeReferencia/SequenciadeReferencia.fa \
     --sjdbGTFfile Dados/SequenciadeReferencia/SequenciadeReferencia.gtf \
     --genomeSAindexNbases 12

# Passo 12: Permitir a execução do script criado no GNU.
chmod +x generate_DadosdoStar.sh

# Passo 13: Instalar o STAR.
wget https://github.com/alexdobin/STAR/archive/refs/tags/2.7.10a.tar.gz

# Passo 14: Extrair o arquivo do STAR.
tar -xzf 2.7.10a.tar.gz

# Passo 15: Navegar até o diretório do STAR.
cd STAR-2.7.10a/bin/Linux_x86_64

# Passo 16: Adicione o STAR ao PAtH
export PATH=$PATH:$(pwd)

# Passo 17: Tornar a mudança permanente
source ~/.bashrc

# Passo 18: Criar ídice de genomas para poder mapear reads com o STAR.
./generate_DadosdoStar.sh

# Passo 19: Criar diretório "MapeamentosFeitos" para salvar os arquivos resultantes do mapeamento do STAR. 
mkdir MapeamentosFeitos

# Passo 20: Abrir o GNU para fazer um script para mapear sequencias de RNA com base na sequencia de referência, utilizando STAR.
nano generate_MapearRNAStar.sh

# Passo 21: Copiar o código a seguir no GNU para poder criar um script de execuçãodo STAR. Para salvar, clicar ctrl+o depois y para aceitar e ctrl+x para fechar o GNU.
#!/bin/bash

STAR --runThreadN 4 \
--readFilesIn Dados/ArquivosParaAnalise/Dados1.fastqsanger Dados/ArquivosParaAnalise/Dados2.fastqsanger \
--genomeDir Dados/DadosdoStar \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix MapeamentosFeitos/Dados1e2

# Passo 22: Permitir a execução do script criado no GNU.
chmod +x generate_MapearRNAStar.sh

# Passo 22: Rodar o script.
./generate_MapearRNAStar.sh

# Passo 23: Utilizar o Multiqc para fazer gráficos a partir de um arquivo gerado pelo STAR. Se der certo, haverá um arquivo em html com os gráficos.
multiqc MapeamentosFeitos/ARQUIVO COM FINAL LOG.FINAL.OUT

# Passo 24: Inspecionar o arquivo em BAM genado pelo STAR. Primeiro é preciso indexar o arquivo utilizando uma ferramenta chamada Samtools. Intalar o Samtools.
sudo apt install samtools

# Passo 25: O arquivo em BAM estará no diretório "MapeamentosFeitos", mas não aparece na interface gráfica até rodar o samtools, então tem que ser pelo terminal mesmo. Basta dar um “ls MapeamentosFeitos” que ele deve aparecer na lista.
samtools index MapeamentosFeitos/ARQUIVO COM FINAL OUT.BAM

# Passo 25: Instalar um programa chamado IGV para seguir com a inspeção
sudo apt install igv

# Passo 26: Abrir o programa IGV para seguir com a inspeção
igv

# Passo 27: Com o IGV aberto, primeiro é preciso carregar a sequencia de referência em FASTA [genomes>load genomes from file>Dados>SequenciadeReferencia>ARQUIVO EM FASTA]. Depois é preciso adicionar os arquivos em BAM [file>load from file>MapeamentosFeitos>ARQUIVOS EM BAM]
Pronto!
