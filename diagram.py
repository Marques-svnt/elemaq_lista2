import matplotlib.pyplot as plt
import numpy as np # Adicionado para manipulação, se necessário

def ler_dados_arquivo_eixo(nome_arquivo):
    """
    Lê o arquivo de dados do eixo gerado pelo C.
    Formato esperado no arquivo:
    # Comentários ou cabeçalho ignorados
    Posicao(m) Torque(N-m) Momento_Fletor(N-m) Forca_Cortante(N)
    """
    dados = {
        "posicao": [],
        "torque": [],
        "momento_fletor": [],
        "forca_cortante": []
    }
    try:
        with open(nome_arquivo, 'r') as f:
            linhas = f.readlines()

        for linha_idx, linha in enumerate(linhas):
            if linha.startswith("#") or linha.strip() == "":
                continue  # ignora cabeçalho e linhas em branco
            
            partes = linha.split()
            if len(partes) == 4:
                try:
                    dados["posicao"].append(float(partes[0]))
                    dados["torque"].append(float(partes[1]))
                    dados["momento_fletor"].append(float(partes[2]))
                    dados["forca_cortante"].append(float(partes[3]))
                except ValueError:
                    print(f"Aviso: Ignorando linha mal formatada no arquivo (linha {linha_idx + 1}): {linha.strip()}")
            else:
                print(f"Aviso: Ignorando linha com número incorreto de colunas (linha {linha_idx + 1}): {linha.strip()}")
                
    except FileNotFoundError:
        print(f"ERRO: Arquivo '{nome_arquivo}' não encontrado.")
        return None
    except Exception as e:
        print(f"ERRO ao ler o arquivo '{nome_arquivo}': {e}")
        return None

    # Para que os diagramas de "degrau" funcionem bem com plt.plot,
    # pode ser necessário duplicar pontos em X para criar as linhas verticais.
    # Ex: se cortante muda de V1 para V2 em X1, precisamos de (X1,V1) e (X1,V2).
    # A lógica de geração em C deve idealmente cuidar disso.
    # Se não, algum pós-processamento aqui pode ser necessário dependendo da saída do C.

    return dados

def plotar_diagramas_eixo(dados, comprimento_eixo=None):
    """
    Plota os diagramas de esforços do eixo.
    """
    if dados is None or not dados["posicao"]:
        print("Não há dados para plotar.")
        return

    x_pos = dados["posicao"]
    torque = dados["torque"]
    momento_fletor = dados["momento_fletor"]
    forca_cortante = dados["forca_cortante"]

    # Determinar o comprimento do eixo para o limite do gráfico se não fornecido
    if comprimento_eixo is None:
        comprimento_eixo = x_pos[-1] if x_pos else 1.0


    fig, axs = plt.subplots(3, 1, figsize=(12, 10), sharex=True)

    # Diagrama de Força Cortante
    axs[0].plot(x_pos, forca_cortante, 'b-', marker='.', linestyle='-', drawstyle='steps-post', label='Força Cortante')
    axs[0].fill_between(x_pos, forca_cortante, step='post', color='blue', alpha=0.2)
    axs[0].set_ylabel("Força Cortante (N)")
    axs[0].set_title("Diagrama de Força Cortante (DFC)")
    axs[0].grid(True, linestyle='--', alpha=0.7)
    axs[0].axhline(0, color='black', linewidth=0.8)

    # Diagrama de Momento Fletor
    axs[1].plot(x_pos, momento_fletor, 'r-', marker='.', linestyle='-', label='Momento Fletor')
    axs[1].fill_between(x_pos, momento_fletor, color='red', alpha=0.2)
    axs[1].set_ylabel("Momento Fletor (N·m)")
    axs[1].set_title("Diagrama de Momento Fletor (DMF)")
    axs[1].grid(True, linestyle='--', alpha=0.7)
    axs[1].axhline(0, color='black', linewidth=0.8)

    # Diagrama de Momento Torsor
    axs[2].plot(x_pos, torque, 'g-', marker='.', linestyle='-', drawstyle='steps-post', label='Momento Torsor')
    axs[2].fill_between(x_pos, torque, step='post', color='green', alpha=0.2)
    axs[2].set_ylabel("Momento Torsor (N·m)")
    axs[2].set_title("Diagrama de Momento Torsor (DMT)")
    axs[2].grid(True, linestyle='--', alpha=0.7)
    axs[2].axhline(0, color='black', linewidth=0.8)

    plt.xlabel("Posição ao longo do eixo (m)")
    
    for ax in axs:
        ax.set_xlim(0, comprimento_eixo)
        # ax.legend() # Opcional

    plt.tight_layout()
    plt.show()

# --- Execução Principal do Script Python ---
if __name__ == "__main__":
    nome_do_arquivo_entrada = "dados_diagramas.txt"
    
    # Tenta ler o comprimento do eixo do cabeçalho do arquivo, se existir
    comprimento_eixo_lido = None
    try:
        with open(nome_do_arquivo_entrada, 'r') as f_temp:
            for linha_temp in f_temp:
                if "# Dados dos Diagramas de Esforços do Eixo (Comprimento:" in linha_temp:
                    comprimento_eixo_lido = float(linha_temp.split("(Comprimento:")[1].split("m)")[0].strip())
                    break
    except:
        pass # Ignora erros ao tentar ler o comprimento

    dados_lidos = ler_dados_arquivo_eixo("dados_diagramas.txt")
    
    if dados_lidos:
        print(f"Dados lidos do arquivo '{nome_do_arquivo_entrada}'. Plotando diagramas...")
        plotar_diagramas_eixo(dados_lidos, comprimento_eixo=comprimento_eixo_lido)
    else:
        print(f"Não foi possível ler ou processar o arquivo '{nome_do_arquivo_entrada}'.")
        print("Certifique-se que o programa C foi executado e gerou o arquivo corretamente.")
        print("Formato esperado no arquivo (após comentários '#'):")
        print("Posicao(m) Torque(N-m) Momento_Fletor(N-m) Forca_Cortante(N)")

