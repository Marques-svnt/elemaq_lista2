# Programa de Cálculo para Análise de Eixos Mecânicos

**UNIVERSIDADE ESTADUAL DE SANTA CRUZ - UESC**
**DEPARTAMENTO DE ENGENHARIAS E COMPUTAÇÃO - DEC**
**CURSO DE ENGENHARIA MECÂNICA**

**Disciplina:** Elementos de Máquina
**Código:** CET 948
**Professor:** Dr. José Carlos de Camargo

---

Este programa realiza a análise de diferentes configurações de eixos mecânicos, conforme especificado na Lista de Exercício 02 - "Projeto de eixos e chavetas". Ele calcula torques, forças nos elementos de transmissão, reações nos mancais e gera dados para os diagramas de esforços solicitantes.

## Pré-requisitos

1.  **Compilador C:** Um compilador C (como GCC/MinGW) é necessário para compilar os arquivos fonte (`main.c`, `q1.c`, `q2.c`, `q3.c`, `q4.c`) e gerar o executável `programa.exe`.
2.  **Python 3:** Necessário para executar o script de plotagem dos diagramas.
3.  **Bibliotecas Python:** As bibliotecas `matplotlib` e `numpy` são usadas para a geração dos gráficos. Para instalá-las, execute o seguinte comando no seu terminal (dentro do diretório do projeto ou onde o `requirements.txt` estiver localizado, se você o colocar lá):
    ```bash
    pip install -r requirements.txt
    ```
    O arquivo `requirements.txt` deve conter:
    ```text
    matplotlib
    numpy
    ```

## Instruções de Execução

1.  **Compilação do Código C:**
    Executar setup.bat

2.  **Executando o Programa:**
    * Após a compilação bem-sucedida, execute `programa.exe`.
    * No Windows, você pode clicar duas vezes no arquivo ou executá-lo via Prompt de Comando:
        ```bash
        programa.exe
        ```
    * Em sistemas Linux/macOS (se compilado para esses sistemas):
        ```bash
        ./programa.exe
        ```

3.  **Saídas do Programa:**
    * O `programa.exe` criará (ou sobrescreverá) a pasta `output/` e, dentro dela:
        * `respostas.txt`: Contém as respostas textuais detalhadas, cálculos e resultados para cada questão analisada.
        * `dados.txt`: Contém os dados numéricos formatados para cada questão, utilizados pelo script Python para gerar os diagramas.
    * Em seguida, o `programa.exe` chamará o script `python src/diagram.py output/dados.txt`.
    * O script Python criará (ou sobrescreverá) a pasta `diagramas/` e salvará as imagens dos diagramas de esforços (ex: `Questao_1_Diagramas_Esforcos.png`, etc.) nesta pasta.

## Funcionamento Detalhado

1.  O `programa.exe` é iniciado.
2.  Abre os arquivos `respostas.txt` e `dados.txt` para escrita.
3.  Processa cada uma das questões (`q1` a `q4`):
    * Calcula torques, forças, reações e outros dados relevantes.
    * Escreve os resultados textuais em `output/respostas.txt`.
    * Escreve os dados numéricos para os diagramas em `output/dados.txt`, formatados conforme esperado pelo script Python.
4.  Fecha os arquivos `respostas.txt` e `dados.txt`.
5.  Chama o script `src/diagram.py`, passando o caminho `output/dados.txt` como argumento.
6.  O script `diagram.py`:
    * Lê os dados de `output/dados.txt`.
    * Cria o diretório `diagramas/` (se não existir).
    * Para cada questão encontrada no arquivo de dados, gera os 5 diagramas de esforços (Cortante Vy, Momento Mz, Cortante Vz, Momento My, Torque T).
    * Salva cada conjunto de diagramas como uma imagem PNG na pasta `diagramas/`.

Certifique-se de que o script `src/diagram.py` está configurado corretamente para ler o arquivo de dados do caminho fornecido como argumento e para salvar as imagens na pasta `diagramas/`.
