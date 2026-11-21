"""
ADAPTADOR: o UNICO arquivo que conhece a estrutura interna do FASST em Python.

A FMU (greenroof_fmu.py) so chama:
    a = FasstAdaptador(**parametros)
    saidas = a.passo(dt, **entradas)

Tarefa do Lucas e do Allan: reescrever os dois metodos abaixo chamando o
FASST transcrito, em vez do ModeloSubstituto.
  - __init__ : o que o programa principal do FASST faz ANTES do laco no tempo
               (ler parametros, montar camadas, condicao inicial).
               NAO ler arquivo meteorologico aqui.
  - passo    : o CORPO do laco no tempo, recebendo a meteorologia como argumento
               e usando T_laje (vinda do DOMUS) como condicao de contorno inferior.

Convencoes do FASST que o adaptador precisa respeitar (conferidas no codigo Fortran):
  1. CONTORNO INFERIOR. O FASST original usa fluxo de calor FIXO de 0,15 W/m2 no fundo
     e o lencol freatico como contorno hidrico. Para o TESTE DE REGRESSAO contra o FASST,
     use contorno="fasst_original" (reproduz isso). Para o acoplamento, contorno="laje".
  2. METEOROLOGIA DO FIM DO PASSO. O FASST avanca [t, t+dt] com a meteorologia de t+dt.
     O mestre FMI precisa entregar em doStep(t) os valores de t+dt (o DOMUS deve setar as
     entradas do fim do passo; nos testes com fmpy a serie e deslocada de um passo).
  3. CHUVA. O FASST le lamina por passo [mm/passo]. A FMU recebe taxa [kg/(m2 s)];
     converter com  mm_no_passo = chuva * dt.
  4. PASSO INTERNO. Piso de 300 s quando o passo meteorologico >= 1 h; no maximo 10 s abaixo
     disso. Com passo de comunicacao de 1-5 min o custo sobe (ver guia, etapa 6).
"""
from fasst.modelo_substituto import ModeloSubstituto


class FasstAdaptador:
    def __init__(self, LAI, h_sub, theta0, T0, contorno="laje"):
        self.contorno = contorno            # "laje" (acoplamento) ou "fasst_original" (regressao)
        # >>> TROCAR pelo init do FASST
        self.m = ModeloSubstituto(LAI=LAI, h_sub=h_sub, theta0=theta0, T0=T0)

    def passo(self, dt, T_ar, UR, G_sol, L_ceu, vento, chuva, T_laje):
        # >>> TROCAR pelo corpo do laco do FASST
        # chuva_mm_passo = chuva * dt     # o FASST espera mm por passo
        # if self.contorno == "fasst_original": fundo com q = 0,15 W/m2 e lencol freatico
        # else: fundo acoplado a laje (T_laje) + membrana/drenagem
        # Deve devolver: q_laje [W/m2, >0 entrando na laje], T_folha [K], T_sub [K],
        #                theta [m3/m3], LE [W/m2], escoamento [kg/m2s]
        return self.m.passo(dt, T_ar, UR, G_sol, L_ceu, vento, chuva, T_laje)
