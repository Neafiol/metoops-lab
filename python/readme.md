
 Laba 1
----

0. Перейдите в папку `python`
1. Для перекомпиляции кода на не линуксе выполните
`gcc -shared -Wl,-install_name,mo.co -o mo.co -fPIC MO2.c`
2. Скачайте зависимости `pip install -r requerements.txt `
3. Запустите GUI `python3 gui.py`

Файлы к первой лабе
- `grafics.py` - рисование графиков для отчета
- `gui.py` - ГУИ
- `MO2.co` - скомпилированные функции под С
- `MO2.c` -исходный код функций переписанный на C

