import matplotlib.pyplot as plt
import numpy as np
import os
import sys # Para ler argumentos da linha de comando

# Funcao ler_dados_eixo_completo (COPIE A ÚLTIMA VERSÃO QUE FUNCIONAVA BEM PARA LEITURA)
def ler_dados_eixo_completo(nome_arquivo):
    dados = {
        "posicao": [], "torque": [],
        "fletor_mz": [], "cortante_vy": [],
        "fletor_my": [], "cortante_vx": [] 
    }
    questoes_lidas = {}
    dados_questao_atual = dados.copy()
    num_questao_atual = 0 

    try:
        with open(nome_arquivo, 'r', encoding='utf-8') as f:
            linhas = f.readlines()

        for linha_idx, linha in enumerate(linhas):
            linha_strip = linha.strip()
            if not linha_strip: continue

            if linha_strip.startswith("# Questao"):
                if num_questao_atual > 0 and dados_questao_atual["posicao"]:
                    questoes_lidas[num_questao_atual] = dados_questao_atual
                elif num_questao_atual == 0 and not questoes_lidas and dados_questao_atual["posicao"]:
                     questoes_lidas[1] = dados_questao_atual
                
                try:
                    num_str = linha_strip.split("Questao")[1].strip().split()[0]
                    num_questao_atual = int(num_str)
                except (IndexError, ValueError) as e:
                    print(f"Aviso (ler_dados): Não foi possível extrair o número da questão da linha: '{linha_strip}'. Erro: {e}")
                    num_questao_atual = -1 
                
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
                    dados_questao_atual["cortante_vx"].append(float(partes[5]))
                except ValueError:
                    print(f"Aviso (ler_dados): Ignorando linha mal formatada (linha {linha_idx + 1}): {linha_strip}")
            elif partes:
                print(f"Aviso (ler_dados): Linha com número incorreto de colunas (linha {linha_idx + 1}): {linha_strip}")
        
        if num_questao_atual != -1 and dados_questao_atual["posicao"]:
            questoes_lidas[num_questao_atual] = dados_questao_atual
        elif not questoes_lidas and dados_questao_atual["posicao"]: 
            questoes_lidas[1] = dados_questao_atual

    except FileNotFoundError:
        print(f"ERRO (ler_dados): Arquivo '{nome_arquivo}' não encontrado.")
        return None
    except Exception as e:
        print(f"ERRO (ler_dados) ao ler o arquivo '{nome_arquivo}': {e}")
        return None
    
    if not questoes_lidas:
        print(f"Nenhum dado de questão válido foi lido do arquivo '{nome_arquivo}'.")
        return None
    return questoes_lidas

# Funcao plotar_diagramas_por_questao (COPIE A ÚLTIMA VERSÃO COM UNIDADES SI E AJUSTES DE FIGSIZE/XLIM)
def plotar_diagramas_por_questao(num_questao, dados, output_dir, comprimento_eixo_calculado=None):
    if not dados["posicao"]:
        print(f"Questão {num_questao}: Sem dados para plotar.")
        return

    x_pos = np.array(dados["posicao"])
    
    unidade_forca = "N"
    unidade_momento = "N·m"
    unidade_posicao = "m"

    fig, axs = plt.subplots(5, 1, figsize=(10, 14), sharex=True) 
    fig.suptitle(f"Diagramas de Esforços (SI) - Questão {num_questao}", fontsize=16, y=0.985)

    axs[0].plot(x_pos, dados["cortante_vy"], 'b-', marker='.', markersize=5, drawstyle='steps-post', label='Vy')
    axs[0].fill_between(x_pos, dados["cortante_vy"], step='post', color='blue', alpha=0.2)
    axs[0].set_ylabel(f"Cortante Vy ({unidade_forca})")
    axs[0].set_title("Força Cortante - Plano Vertical (DFC-y)")
    axs[0].grid(True, linestyle=':', alpha=0.6); axs[0].axhline(0, color='black', lw=0.7)

    axs[1].plot(x_pos, dados["fletor_mz"], 'r-', marker='.', markersize=5, label='Mz')
    axs[1].fill_between(x_pos, dados["fletor_mz"], color='red', alpha=0.2)
    axs[1].set_ylabel(f"Momento Fletor Mz ({unidade_momento})")
    axs[1].set_title("Momento Fletor - Plano Vertical (DMF-z)")
    axs[1].grid(True, linestyle=':', alpha=0.6); axs[1].axhline(0, color='black', lw=0.7)

    axs[2].plot(x_pos, dados["cortante_vx"], 'c-', marker='.', markersize=5, drawstyle='steps-post', label='Vz (Vx no arquivo)')
    axs[2].fill_between(x_pos, dados["cortante_vx"], step='post', color='cyan', alpha=0.2)
    axs[2].set_ylabel(f"Cortante Vz ({unidade_forca})")
    axs[2].set_title("Força Cortante - Plano Horizontal (DFC-z)")
    axs[2].grid(True, linestyle=':', alpha=0.6); axs[2].axhline(0, color='black', lw=0.7)

    axs[3].plot(x_pos, dados["fletor_my"], 'm-', marker='.', markersize=5, label='My')
    axs[3].fill_between(x_pos, dados["fletor_my"], color='magenta', alpha=0.2)
    axs[3].set_ylabel(f"Momento Fletor My ({unidade_momento})")
    axs[3].set_title("Momento Fletor - Plano Horizontal (DMF-y)")
    axs[3].grid(True, linestyle=':', alpha=0.6); axs[3].axhline(0, color='black', lw=0.7)

    axs[4].plot(x_pos, dados["torque"], 'g-', marker='.', markersize=5, drawstyle='steps-post', label='T')
    axs[4].fill_between(x_pos, dados["torque"], step='post', color='green', alpha=0.2)
    axs[4].set_ylabel(f"Momento Torsor T ({unidade_momento})")
    axs[4].set_title("Diagrama de Momento Torsor (DMT)")
    axs[4].grid(True, linestyle=':', alpha=0.6); axs[4].axhline(0, color='black', lw=0.7)

    plt.xlabel(f"Posição ao longo do eixo ({unidade_posicao})", fontsize=11)
    
    if len(x_pos) > 0:
        min_data_x = min(x_pos)
        max_data_x = max(x_pos)
        if comprimento_eixo_calculado is not None and comprimento_eixo_calculado > max_data_x:
            max_data_x = comprimento_eixo_calculado
        plot_range_x = max_data_x - min_data_x
        padding_x = plot_range_x * 0.05 if plot_range_x > 1e-6 else 0.1 
        lim_min_x = min_data_x - padding_x
        lim_max_x = max_data_x + padding_x
        if min_data_x >= 0 and lim_min_x < 0:
            lim_min_x = -padding_x 
            if max_data_x == 0: lim_max_x = padding_x
    else: 
        lim_min_x = 0; lim_max_x = 1 
            
    for ax in axs:
        ax.set_xlim(lim_min_x, lim_max_x) 
        ax.tick_params(axis='y', labelsize=8)
        ax.tick_params(axis='x', labelsize=8)
        ax.title.set_fontsize(10)
        ax.yaxis.label.set_fontsize(9)

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])

    nome_arquivo_plot = f"Questao_{num_questao}_Diagramas_SI.png"
    caminho_completo = os.path.join(output_dir, nome_arquivo_plot)
    try:
        plt.savefig(caminho_completo, dpi=150)
        print(f"Diagrama da Questão {num_questao} (SI) salvo em: {caminho_completo}")
    except Exception as e:
        print(f"ERRO ao salvar o diagrama da Questão {num_questao} em '{caminho_completo}': {e}")
    
    plt.close(fig)

if __name__ == "__main__":
    # Prioriza o argumento de linha de comando (passado pelo C)
    if len(sys.argv) > 1:
        nome_arquivo_dados = sys.argv[1]
        print(f"Lendo dados de: '{nome_arquivo_dados}' (fornecido como argumento)")
    else:
        # Fallback se o script for executado manualmente sem argumentos
        # Assume que o script está em 'src/' e 'output/' está no mesmo nível do diretório pai de 'src/'
        # Ou, se executado da raiz do projeto, 'output/dados.txt'
        # Para maior robustez, o C deve SEMPRE passar o argumento.
        # Este fallback assume que o script é executado da raiz do projeto.
        nome_arquivo_dados = os.path.join("output", "dados.txt")
        print(f"AVISO: Nenhum caminho de arquivo de dados fornecido como argumento.")
        print(f"       Tentando o caminho padrão: '{nome_arquivo_dados}'")
        # Verificação adicional se o caminho padrão existe, para ajudar na depuração
        if not os.path.exists(nome_arquivo_dados):
             alt_path = "dados.txt" # Tenta no diretório atual do script
             print(f"AVISO: Caminho padrão '{nome_arquivo_dados}' não encontrado. Tentando '{alt_path}' no diretório atual.")
             if os.path.exists(alt_path):
                 nome_arquivo_dados = alt_path
             else:
                 print(f"ERRO: Arquivo de dados não encontrado em '{nome_arquivo_dados}' nem em '{alt_path}'.")
                 sys.exit(1) # Sai se não encontrar o arquivo de dados no fallback

    diretorio_saida_diagramas = "diagramas"

    if not os.path.exists(diretorio_saida_diagramas):
        try:
            os.makedirs(diretorio_saida_diagramas)
            print(f"Diretório '{diretorio_saida_diagramas}' criado com sucesso.")
        except OSError as e:
            print(f"ERRO ao criar o diretório '{diretorio_saida_diagramas}': {e}")
            diretorio_saida_diagramas = "." # Fallback

    dados_por_questao = ler_dados_eixo_completo(nome_arquivo_dados)

    if dados_por_questao:
        for num_q, dados_q in dados_por_questao.items():
            if num_q == -1:
                print("Aviso: Ignorando plotagem para questão com número indefinido.")
                continue
            print(f"\nProcessando e plotando diagramas (SI) para Questão {num_q}...")
            
            compr_eixo_final = None
            if dados_q["posicao"]:
                min_pos_dados = min(dados_q["posicao"])
                max_pos_dados = max(dados_q["posicao"])
                range_dados = max_pos_dados - min_pos_dados
                
                # Define o comprimento do eixo para plotagem para cobrir o intervalo dos dados
                # com uma pequena margem, ou um mínimo se o range for muito pequeno.
                if range_dados < 1e-6 : # Se todos os pontos x são (quase) os mesmos
                    compr_eixo_final = max_pos_dados + 0.1 # Adiciona uma pequena extensão
                else:
                    compr_eixo_final = max_pos_dados + range_dados * 0.05 # 5% de padding no final
                
                # Garante que o comprimento do eixo seja pelo menos um valor mínimo razoável
                compr_eixo_final = max(compr_eixo_final, 0.1) 


            plotar_diagramas_por_questao(num_q, dados_q, diretorio_saida_diagramas, compr_eixo_final)
    else:
        print(f"Não foi possível ler ou processar dados do arquivo '{nome_arquivo_dados}'.")