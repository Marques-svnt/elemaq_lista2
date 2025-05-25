@echo off
echo Verificando e instalando dependencias Python...
pip install -r src/requirements.txt
if %errorlevel% neq 0 (
    echo ERRO: Falha ao instalar dependencias Python.
    echo Certifique-se de que Python e pip estao instalados e no PATH.
    pause
    exit /b 1
)
echo Dependencias Python instaladas/verificadas com sucesso.
pause