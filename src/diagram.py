import matplotlib.pyplot as plt
import numpy as np
import os # Importa o módulo os

def ler_dados_eixo_completo(nome_arquivo):
    dados = {
        "posicao": [], "torque": [],
        "fletor_mz": [], "cortante_vy": [],
        "fletor_my": [], "cortante_vx": [] # Lembre-se que vx aqui corresponde a Vz
    }
    questoes_lidas = {}
    dados_questao_atual = dados.copy()
    num_questao_atual = 0

    try:
        with open(nome_arquivo, 'r') as f:
            linhas = f.readlines()

        for linha_idx, linha in enumerate(linhas):
            linha_strip = linha.strip()
            if not linha_strip: continue

            if linha_strip.startswith("# Questao"):
                if num_questao_atual > 0 and dados_questao_atual["posicao"]:
                    questoes_lidas[num_questao_atual] = dados_questao_atual
                
                try:
                    # Extrai o número da questão de forma mais robusta
                    num_str = linha_strip.split("Questao")[1].strip().split()[0]
                    num_questao_atual = int(num_str)
                except (IndexError, ValueError) as e:
                    print(f"Aviso: Não foi possível extrair o número da questão da linha: '{linha_strip}'. Erro: {e}")
                    num_questao_atual = -1 # Marca como indefinida
                
                dados_questao_atual = {
                    "posicao": [], "torque": [],
                    "fletor_mz": [], "cortante_vy": [],
                    "fletor_my": [], "cortante_vx": []
                }
                continue
            
            if linha_strip.startswith("#"): continue

            partes = linha_strip.split()
            if len(partes) == 6:
                try:
                    dados_questao_atual["posicao"].append(float(partes[0]))
                    dados_questao_atual["torque"].append(float(partes[1]))
                    dados_questao_atual["fletor_mz"].append(float(partes[2]))
                    dados_questao_atual["cortante_vy"].append(float(partes[3]))
                    dados_questao_atual["fletor_my"].append(float(partes[4]))
                    dados_questao_atual["cortante_vx"].append(float(partes[5])) # Vz
                except ValueError:
                    print(f"Aviso: Ignorando linha mal formatada (linha {linha_idx + 1}): {linha_strip}")
            elif partes:
                print(f"Aviso: Linha com número incorreto de colunas (linha {linha_idx + 1}): {linha_strip}")
        
        if num_questao_atual != -1 and dados_questao_atual["posicao"]: # Salva a última questão lida
            questoes_lidas[num_questao_atual] = dados_questao_atual
        elif not questoes_lidas and dados_questao_atual["posicao"]:
            questoes_lidas[1] = dados_questao_atual # Caso sem header de questão

    except FileNotFoundError:
        print(f"ERRO: Arquivo '{nome_arquivo}' não encontrado.")
        return None
    except Exception as e:
        print(f"ERRO ao ler o arquivo '{nome_arquivo}': {e}")
        return None
    
    if not questoes_lidas:
        print("Nenhum dado de questão foi lido do arquivo.")
        return None

    return questoes_lidas

# Modificada para aceitar o diretório de saída
def plotar_diagramas_por_questao(num_questao, dados, output_dir, comprimento_eixo=None):
    if not dados["posicao"]:
        print(f"Questão {num_questao}: Sem dados para plotar.")
        return

    x_pos = np.array(dados["posicao"])
    
    if comprimento_eixo is None and len(x_pos) > 0:
        comprimento_eixo = x_pos[-1]
    elif comprimento_eixo is None:
        comprimento_eixo = 10.0 # Um valor padrão se não houver dados
    
    # Define as unidades com base no problema (lb e in)
    unidade_forca = "lb"
    unidade_momento = "lb·in"
    unidade_posicao = "in"


    fig, axs = plt.subplots(5, 1, figsize=(12, 20), sharex=True) # Aumentei um pouco o figsize
    fig.suptitle(f"Diagramas de Esforços - Questão {num_questao}", fontsize=18, y=0.99)

    # Cortante Vy (Plano Vertical)
    axs[0].plot(x_pos, dados["cortante_vy"], 'b-', marker='.', drawstyle='steps-post', label='Vy')
    axs[0].fill_between(x_pos, dados["cortante_vy"], step='post', color='blue', alpha=0.2)
    axs[0].set_ylabel(f"Cortante Vy ({unidade_forca})")
    axs[0].set_title("Força Cortante - Plano Vertical (DFC-y)")
    axs[0].grid(True); axs[0].axhline(0, color='black', lw=0.5)

    # Momento Mz (Plano Vertical)
    axs[1].plot(x_pos, dados["fletor_mz"], 'r-', marker='.', label='Mz')
    axs[1].fill_between(x_pos, dados["fletor_mz"], color='red', alpha=0.2)
    axs[1].set_ylabel(f"Momento Fletor Mz ({unidade_momento})")
    axs[1].set_title("Momento Fletor - Plano Vertical (DMF-z)")
    axs[1].grid(True); axs[1].axhline(0, color='black', lw=0.5)

    # Cortante Vz (Plano Horizontal - lido como 'cortante_vx')
    axs[2].plot(x_pos, dados["cortante_vx"], 'c-', marker='.', drawstyle='steps-post', label='Vz (Vx no arquivo)')
    axs[2].fill_between(x_pos, dados["cortante_vx"], step='post', color='cyan', alpha=0.2)
    axs[2].set_ylabel(f"Cortante Vz ({unidade_forca})")
    axs[2].set_title("Força Cortante - Plano Horizontal (DFC-z)")
    axs[2].grid(True); axs[2].axhline(0, color='black', lw=0.5)

    # Momento My (Plano Horizontal)
    axs[3].plot(x_pos, dados["fletor_my"], 'm-', marker='.', label='My')
    axs[3].fill_between(x_pos, dados["fletor_my"], color='magenta', alpha=0.2)
    axs[3].set_ylabel(f"Momento Fletor My ({unidade_momento})")
    axs[3].set_title("Momento Fletor - Plano Horizontal (DMF-y)")
    axs[3].grid(True); axs[3].axhline(0, color='black', lw=0.5)

    # Momento Torsor T
    axs[4].plot(x_pos, dados["torque"], 'g-', marker='.', drawstyle='steps-post', label='T')
    axs[4].fill_between(x_pos, dados["torque"], step='post', color='green', alpha=0.2)
    axs[4].set_ylabel(f"Momento Torsor T ({unidade_momento})")
    axs[4].set_title("Diagrama de Momento Torsor (DMT)")
    axs[4].grid(True); axs[4].axhline(0, color='black', lw=0.5)

    plt.xlabel(f"Posição ao longo do eixo ({unidade_posicao})")
    if comprimento_eixo is not None:
        for ax in axs:
            ax.set_xlim(0, comprimento_eixo) # Garante que o eixo x comece em 0
            
    plt.tight_layout(rect=[0, 0.03, 1, 0.97]) # Ajusta para o supertítulo e eixos

    # Salva a figura
    nome_arquivo_plot = f"questao_{num_questao}_diagramas.png"
    caminho_completo = os.path.join(output_dir, nome_arquivo_plot)
    try:
        plt.savefig(caminho_completo, dpi=150) # dpi aumenta a resolução da imagem salva
        print(f"Diagrama da Questão {num_questao} salvo em: {caminho_completo}")
    except Exception as e:
        print(f"ERRO ao salvar o diagrama da Questão {num_questao}: {e}")
    
    plt.close(fig) # Fecha a figura para liberar memória

if __name__ == "__main__":
    nome_arquivo_dados = "dados.txt"
    diretorio_saida_diagramas = "diagramas" # Define o nome da pasta

    # Cria o diretório de saída se não existir
    if not os.path.exists(diretorio_saida_diagramas):
        try:
            os.makedirs(diretorio_saida_diagramas)
            print(f"Diretório '{diretorio_saida_diagramas}' criado com sucesso.")
        except OSError as e:
            print(f"ERRO ao criar o diretório '{diretorio_saida_diagramas}': {e}")
            # Decide se quer sair ou tentar salvar no diretório atual
            diretorio_saida_diagramas = "." # Tenta salvar no diretório atual como fallback

    dados_por_questao = ler_dados_eixo_completo(nome_arquivo_dados)

    if dados_por_questao:
        for num_q, dados_q in dados_por_questao.items():
            if num_q == -1 : # Questão com parsing de número falho
                print("Aviso: Ignorando plotagem para questão com número indefinido.")
                continue
            print(f"\nPlotando e salvando diagramas para Questão {num_q}...")
            compr_eixo = None
            if dados_q["posicao"]:
                # Garante que o comprimento do eixo seja pelo menos o máximo x, ou um valor padrão se vazio
                max_pos = max(dados_q["posicao"]) if dados_q["posicao"] else 0.0
                compr_eixo = max(max_pos, 10.0) # Usa um comprimento mínimo para o eixo X do gráfico
            
            plotar_diagramas_por_questao(num_q, dados_q, diretorio_saida_diagramas, compr_eixo)
    else:
        print(f"Não foi possível ler ou processar o arquivo '{nome_arquivo_dados}'.")