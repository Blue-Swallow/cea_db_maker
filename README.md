# RockCombustSim
Rocket Combustion Simulator using  
NASA-CEAを使った簡易ロケット燃焼シミュレータ(仮) ← 本体は未実装  
※ windows, python 3 に対応  
※ 必要なpythonパッケージは適宜入れてください（不親切ですいません）

## How to use  (for busy man)　忙しい人のための使い方解説．  
本ツールを使ってNASA-CEAで計算できるパラメタのデータリストを作成する．  
具体的には，計算したいO/F,Pcの範囲を指定すれば，パラメタ毎に結果がまとまったcsvファイル群を出力してくれる．

### 1. 準備
___
* 適当な場所に本プロジェクトをクローン or ダウンロード
* "[CEA+Fortran.zip](https://www.grc.nasa.gov/WWW/CEAWeb/CEA+Fortran.zip)" と "[CEAexec-win.zip](https://www.grc.nasa.gov/WWW/CEAWeb/CEAexec-win.zip)" を [NASA Chmmical Equilibrium with Application](https://www.grc.nasa.gov/WWW/CEAWeb/ceaguiDownload-win.htm) からダウンロードし展開
* プロジェクト中の".replace"拡張子ファイルを，NASAからダウンロードした同名ファイルと入れ替え，同一フォルダにまとめる

**必要ファイル一覧**   
* execute.py
* gen_inp.py
* FCEA2.exe
* b1b2b3.exe
* syntax.exe
* thermo.inp
* thermo.lib
* trans.inp
* trans.lib

### 2. inp (CEA入力ファイル) の生成
___
* `gen_inp.py` を実行

* 言語の選択　jp=日本語，en=英語
~~~
Please select language.
['jp', 'en']

>> jp
~~~

* 平衡計算条件の指定．（０ or 1 or 2 を入力）
~~~
計算オプション(0~2)を選択してください．
例: 0: 全域平衡計算
    1: 燃焼器内のみ平衡計算
    2: スロートまで平衡計算

>> 0
~~~

* 酸化剤の種類を指定
~~~
酸化剤の種類を入力してください.
*記号はNASA RP-1311-P2 app.B に準拠
例: O2(L)

>> O2(L)
~~~

* 燃料の種類を指定 (半角英字なら何でも可)
~~~
燃料の名前を入力してください
例: PMMA

>> PMMA
~~~

* 酸化剤の初期温度を入力
~~~
酸化剤の初期温度[K]を入力してください

>> 90
~~~

* 燃料の初期温度を入力
~~~
燃料の初期温度[K]を入力してください

>> 280
~~~

* 燃料の標準生成エンタルピを入力
~~~
燃料の標準生成エンタルピ[kJ/mol]を入力してください

>> -468.3
~~~

* 1 mol 燃料中の元素とそのmol数を入力　※各文字間には半角スペースをいれる
~~~
1molの燃料に含まれる元素とそのmol数(整数)を入力してください.
例: C 5 H 2 O 6

>> C 5 H 2 O 6
~~~

* ノズル開口比を入れる　※わからなければ1.0
~~~
開口比Ae/Atを入力してください.

>> 2.0
~~~

* 計算したいO/F の範囲を入力　※各数値間には半角スペースをいれる
~~~
計算するO/Fの範囲を入力してください.
例) 0.5~10 を 0.1毎に計算する場合.
0.5　10　0.1

>> 0.5 10 0.1
~~~

* 計算したい燃焼室圧力の範囲を入力　※各数値間には半角スペースをいれる
~~~
計算する燃焼室圧力[MPa]の範囲を入力してください.
例) 0.5 MPa ~ 5.0 MPa を 0.1 MPa毎に計算する場合.

>> 0.5　5.0　0.1
~~~

* 計算条件ファイル群出力時の格納フォルダ名を入力
~~~
Input a Case Name (Folder Name) >>

>> test
~~~

* 同じディレクトリ上に(今回は)inpファイル群を格納した `test` フォルダが生成される

### 3. CEAの実行
___
* `execute.py` を実行

* 先程生成したフォルダ名を入力
~~~
Input Folder Name (e.g. "O2+PMMA")>>

>> test
~~~

* とりあえず `n` を入力
~~~
Input Polyberization Number. If you didn't assign it, please input "n"
(Directory structure is like "O2+PMMA/n=100")

>> n
~~~

* 100%になって完了すれば，指定したフォルダ (今回は `test` ) 中のフォルダ `csv_database` に各パラメタ毎のcsvファイル，`out`に各計算条件毎のCEA出力(out)ファイルが生成されているはず  
※csvファルは，縦軸O/F,横軸Pc(燃焼室圧力)
  

  
**ロケットパラメタ一覧**

| Symbol | Parameter |
|:---|:---|
|CSTAR |c* 効率 |
|CF |推力係数 |
|Isp |比推力(最適膨張) [s] |
|Ivac |真空中比推力 [s] |
  
**熱力学的パラメタ一覧**

| Symbol | Parameter |
|:---|:---|
|Cp |定圧比熱 [kJ/kg-K] |
|G |ギブス自由エネルギ [kJ/kg] |
|GAMMAs |比熱比 |
|H |エンタルピ [kJ/kg] |
|M |モル質量 [kg/mol] |
|MACH |マッハ数 |
|P |圧力 [MPa] |
|RHO |密度 [kg/m^3] |
|S |エントロピ [kJ/kg-K] |
|SON |音速 [m/s] |
|T |温度 [K] |
|U |内部エネルギ [ｋJ/kg] |

**熱輸送関係のパラメタ一覧**

| Symbol | Parameter |
|:---|:---|
|VISC |粘性係数 [mP] |
|CONDUCTIVITY |熱伝導率 [mW/cm-K] |
|PLANDTL |プラントル数 [-] |
※frozen 条件で計算した**熱伝導率**と**プラントル数**は反応による生成物の影響を考慮していないため低い値が得られる。正しく計算したい場合はequilibrium条件で計算しなければいけない。


**添字一覧**

| Symbol | Subscription |
|:---|:---|
|_c |燃焼室後方における計算結果 |
|_t |ノズルスロートにおける計算結果 |
|_e |ノズル出口における計算結果 |
