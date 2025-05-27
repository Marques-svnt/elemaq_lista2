import os
import re
from fpdf import FPDF # Lembre-se: pip install fpdf2

try:
    from PIL import Image as PILImage
    HAS_PILLOW = True
except ImportError:
    HAS_PILLOW = False
    print("AVISO: Biblioteca Pillow (PIL) não encontrada. O dimensionamento da imagem pode ser menos preciso.")
    print("       Para melhores resultados, instale Pillow: pip install Pillow")

# Constantes
ARQUIVO_RESPOSTAS = "respostas.txt"
PASTA_DIAGRAMAS = "diagramas"
ARQUIVO_PDF_SAIDA = "Relatorio_Analise_Eixos_Finalizado.pdf" # Novo nome para teste

class PDFRelatorio(FPDF):
    def header(self):
        self.set_font('Arial', 'B', 14)
        title = 'Relatório de Análise de Eixos Mecânicos'
        # Para centralizar, podemos calcular a largura do texto e posicionar
        title_w = self.get_string_width(title)
        self.set_x((self.w - title_w) / 2)
        self.cell(title_w, 10, title, 0, 1, 'C')
        self.ln(5)

    def footer(self):
        self.set_y(-15)
        self.set_font('Arial', 'I', 8)
        self.cell(0, 10, f'Página {self.page_no()}/{{nb}}', 0, 0, 'C')

    def titulo_secao_principal(self, titulo):
        self.set_font('Arial', 'B', 13)
        self.ln(5)
        self.cell(0, 8, titulo, 0, 1, 'L')
        self.ln(3)

    def corpo_texto_questao(self, bloco_texto_questao):
        self.set_font('Courier', '', 10.5) # Fonte para o corpo do texto
        largura_texto = self.w - self.l_margin - self.r_margin # Largura útil

        for linha_original in bloco_texto_questao:
            linha_para_escrever = linha_original.rstrip('\n')

            if linha_para_escrever.startswith("--- Questao") or \
               linha_para_escrever.startswith("--- Fim da Questao"):
                # Se o título da questão já foi adicionado por titulo_secao_principal,
                # podemos pular essas linhas aqui, ou formatá-las de forma diferente.
                # Por agora, vamos pular para evitar duplicidade se titulo_secao_principal já as cobre.
                continue
            
            if "Refazendo a Questao" in linha_para_escrever:
                self.ln(1) # Espaço antes
                self.set_font('Courier', 'BI', 10.5) # Negrito e Itálico
                self.set_x(self.l_margin)
                self.multi_cell(largura_texto, 6, linha_para_escrever.encode('latin-1', 'replace').decode('latin-1'), 0, 'L')
                self.set_font('Courier', '', 10.5) # Volta para a fonte normal
                self.ln(1) # Espaço depois
                continue

            try:
                linha_decodificada = linha_para_escrever.encode('latin-1', 'replace').decode('latin-1')
                self.set_x(self.l_margin) 
                self.multi_cell(largura_texto, 6, linha_decodificada, 0, 'L') # Altura da linha 6
            except self.fpdf.errors.FPDFException as e:
                print(f"AVISO FPDFException em multi_cell: {e} - Linha: '{linha_para_escrever[:50]}...'")
                self.set_text_color(255,0,0)
                self.set_x(self.l_margin)
                self.multi_cell(largura_texto, 6, "[ERRO DE RENDERIZAÇÃO DE TEXTO]", 0, 'L')
                self.set_text_color(0,0,0)
            # Não adicionar self.ln() após cada multi_cell, pois ela já faz isso.
        self.ln(4) # Espaço ao final do bloco de texto

    def adicionar_diagrama_questao(self, caminho_diagrama): # Removido nome_arquivo_diagrama como arg, pego de caminho_diagrama
        nome_arquivo_diagrama = os.path.basename(caminho_diagrama)
        if not os.path.exists(caminho_diagrama):
            self.set_font('Arial', 'I', 10); self.set_text_color(255, 0, 0)
            self.multi_cell(self.w - self.l_margin - self.r_margin, 6, f"ERRO: Diagrama '{caminho_diagrama}' não encontrado.")
            self.set_text_color(0, 0, 0); return

        # Título do Diagrama (já estamos em uma nova página dedicada)
        self.set_font('Arial', 'B', 12) # Título um pouco maior
        self.cell(0, 10, f"Diagramas de Esforços: {nome_arquivo_diagrama.replace('_Diagramas_Esforcos.png', '').replace('_', ' ')}", 0, 1, 'C')
        self.ln(3)

        # Dimensionamento da Imagem para caber na página
        y_atual = self.get_y()
        altura_disponivel_total = self.h - y_atual - self.b_margin - 5 # 5mm de folga inferior
        largura_disponivel_total = self.w - self.l_margin - self.r_margin

        if altura_disponivel_total <= 10: # Se não há quase espaço, não tenta desenhar
            print(f"AVISO: Pouco espaço vertical para o diagrama '{nome_arquivo_diagrama}' na página.")
            return

        img_w_pdf, img_h_pdf = largura_disponivel_total, 0 # Padrão: preenche largura, altura proporcional

        if HAS_PILLOW:
            try:
                with PILImage.open(caminho_diagrama) as img_pil:
                    img_pil_w, img_pil_h = img_pil.size
                    if img_pil_w == 0 or img_pil_h == 0: raise ValueError("Dimensões da imagem são zero.")
                    aspect_ratio = float(img_pil_h) / float(img_pil_w)

                    w_calc_por_largura = largura_disponivel_total
                    h_calc_por_largura = w_calc_por_largura * aspect_ratio

                    h_calc_por_altura = altura_disponivel_total
                    w_calc_por_altura = h_calc_por_altura / aspect_ratio
                    
                    # Escolhe o dimensionamento que melhor usa o espaço mantendo a proporção
                    if h_calc_por_largura <= altura_disponivel_total:
                        # Ajustar pela largura cabe na altura
                        img_w_pdf, img_h_pdf = w_calc_por_largura, h_calc_por_largura
                    else:
                        # Ajustar pela altura é necessário
                        img_w_pdf, img_h_pdf = w_calc_por_altura, h_calc_por_altura
            except Exception as e_pil:
                print(f"AVISO: Erro ao usar Pillow para '{caminho_diagrama}': {e_pil}. Usando dimensionamento FPDF padrão.")
                img_w_pdf = largura_disponivel_total; img_h_pdf = 0
        else: # Sem Pillow
            img_w_pdf = largura_disponivel_total
            img_h_pdf = 0 # FPDF tentará manter a proporção, mas pode transbordar a altura

        img_w_pdf = max(1, img_w_pdf)
        if img_h_pdf != 0: img_h_pdf = max(1, img_h_pdf)

        pos_x_img = (self.w - img_w_pdf) / 2 # Centraliza a imagem
        
        try:
            self.image(caminho_diagrama, x=pos_x_img, y=self.get_y(), w=img_w_pdf, h=img_h_pdf)
        except RuntimeError as e:
            self.set_font('Arial', 'I', 10); self.set_text_color(255,0,0)
            self.multi_cell(self.w - self.l_margin - self.r_margin, 6, f"ERRO ao renderizar imagem '{nome_arquivo_diagrama}': {e}.")
            self.set_text_color(0,0,0)


def extrair_numero_questao_do_nome_arquivo(nome_arquivo):
    match = re.search(r'[Qq]uestao_(\d+)', nome_arquivo)
    if match: return int(match.group(1))
    match_simples = re.search(r'q(\d+)', nome_arquivo, re.IGNORECASE)
    if match_simples: return int(match_simples.group(1))
    return None

def parse_respostas_txt(caminho_arquivo):
    respostas_por_questao = {}
    bloco_atual = []
    num_questao_corrente = None
    try:
        with open(caminho_arquivo, 'r', encoding='utf-8') as f:
            linhas_completas = f.readlines()
        for linha in linhas_completas:
            strip_linha = linha.strip()
            match_inicio = re.match(r"--- Questao (\d+) ---", strip_linha)
            match_fim = re.match(r"--- Fim da Questao (\d+) ---", strip_linha)

            if match_inicio:
                if num_questao_corrente is not None and bloco_atual:
                    respostas_por_questao[num_questao_corrente] = bloco_atual
                num_questao_corrente = int(match_inicio.group(1))
                bloco_atual = [] # Começa um novo bloco DEPOIS do marcador de início
            elif match_fim and num_questao_corrente is not None:
                if int(match_fim.group(1)) == num_questao_corrente:
                    if bloco_atual: # Salva o bloco antes de resetar
                        respostas_por_questao[num_questao_corrente] = bloco_atual
                    bloco_atual = []
                    num_questao_corrente = None # Finalizou este bloco de questão
            elif num_questao_corrente is not None: # Estamos dentro de um bloco de questão
                bloco_atual.append(linha) # Adiciona a linha original (com seu \n)
    
        if num_questao_corrente is not None and bloco_atual: # Pega o último bloco se não houver marcador de fim
            respostas_por_questao[num_questao_corrente] = bloco_atual
        
        if not respostas_por_questao and linhas_completas: 
            print("AVISO: Nenhum marcador '--- Questao X ---' encontrado. Tratando todo o conteúdo como Questão 1.")
            respostas_por_questao[1] = linhas_completas

    except FileNotFoundError:
        print(f"ERRO: Arquivo de respostas '{caminho_arquivo}' não encontrado.")
        return None, True
    except Exception as e:
        print(f"ERRO ao processar '{caminho_arquivo}': {e}")
        return None, True
    
    if not respostas_por_questao: print(f"Nenhum bloco de questão lido de '{caminho_arquivo}'.")
    return respostas_por_questao, False


def gerar_pdf():
    pdf = PDFRelatorio(orientation='P', unit='mm', format='A4')
    pdf.alias_nb_pages()

    respostas_por_questao, erro_critico_respostas = parse_respostas_txt(ARQUIVO_RESPOSTAS)

    if erro_critico_respostas:
        pdf.add_page(); pdf.set_font('Arial', 'B', 12); pdf.set_text_color(255,0,0)
        pdf.multi_cell(0, 10, f"ERRO CRÍTICO ao ler ou processar '{ARQUIVO_RESPOSTAS}'.")
        try: pdf.output(ARQUIVO_PDF_SAIDA, "F"); print(f"PDF de erro '{ARQUIVO_PDF_SAIDA}' gerado.")
        except Exception as e_pdf: print(f"ERRO adicional ao salvar PDF de erro: {e_pdf}")
        return

    mapa_diagramas_por_questao = {}
    if not os.path.isdir(PASTA_DIAGRAMAS):
        print(f"AVISO: Pasta de diagramas '{PASTA_DIAGRAMAS}' não encontrada.")
    else:
        for nome_arquivo in sorted(os.listdir(PASTA_DIAGRAMAS)):
            if nome_arquivo.lower().endswith(('.png', '.jpg', '.jpeg')):
                num_q_diag = extrair_numero_questao_do_nome_arquivo(nome_arquivo)
                if num_q_diag is not None:
                    mapa_diagramas_por_questao[num_q_diag] = os.path.join(PASTA_DIAGRAMAS, nome_arquivo)

    todos_numeros_questoes = sorted(list(set(respostas_por_questao.keys()) | set(mapa_diagramas_por_questao.keys())))

    if not todos_numeros_questoes:
        pdf.add_page(); pdf.titulo_secao_principal("Relatório Vazio")
        pdf.set_font('Arial', '', 11); pdf.multi_cell(pdf.w - pdf.l_margin - pdf.r_margin,10,"Nenhum conteúdo encontrado.")
    else:
        for num_q in todos_numeros_questoes:
            # --- Página para Texto da Questão ---
            pdf.add_page() 
            pdf.titulo_secao_principal(f"Análise da Questão {num_q}")
            if num_q in respostas_por_questao:
                pdf.corpo_texto_questao(respostas_por_questao[num_q])
            else: 
                pdf.set_font('Arial', 'I', 10)
                pdf.multi_cell(pdf.w - pdf.l_margin - pdf.r_margin, 6, "Nenhum texto de resposta detalhada encontrado.")
                pdf.ln(5)

            # --- Página para Diagrama da Questão ---
            if num_q in mapa_diagramas_por_questao:
                pdf.add_page() # <<--- SEMPRE NOVA PÁGINA PARA O DIAGRAMA
                pdf.adicionar_diagrama_questao(mapa_diagramas_por_questao[num_q])
            elif num_q in respostas_por_questao: 
                # Se houve texto mas não diagrama, não precisa de mensagem adicional na página do diagrama.
                pass
            
    try:
        pdf.output(ARQUIVO_PDF_SAIDA, "F")
        print(f"Relatório PDF '{ARQUIVO_PDF_SAIDA}' gerado com sucesso!")
    except Exception as e:
        print(f"ERRO ao salvar o PDF '{ARQUIVO_PDF_SAIDA}': {e}")

if __name__ == '__main__':
    gerar_pdf()