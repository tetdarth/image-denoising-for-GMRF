# ガウス型マルコフ確率場を用いた画像ノイズ除去
## IVGMRF
確率変数間の関係を表すグラフィカルモデルである**マルコフ確率場**を用いて画像ノイズ除去を行います。
特に、各確率変数が多次元ガウス分布に従って生成されると仮定したのがガウス型マルコフ確率場であり、DNNを用いたノイズ除去に比べて高速であり、教師なし学習なので未知のデータに対しても対応することができます。
(3枚の劣化画像から修復)
![ivgmrf](https://github.com/tetdarth/image-denoising-for-GMRF/assets/136053901/0e43169f-75ee-448e-af66-c3bc1b09315c)

## DVGMRF
IVGMRFにおける尤度の設計を拡張したモデルです。
![dvgmrf](https://github.com/tetdarth/image-denoising-for-GMRF/assets/136053901/e5a9e802-70a6-409c-9910-295a28a4fccf)

## IVHGMRF
IVGMRFにおけるバイアスベクトルをガウス分布に従うと仮定した上で拡張したモデルです。
![ivhgmrf](https://github.com/tetdarth/image-denoising-for-GMRF/assets/136053901/438aea57-c813-44fa-9c0b-bcd42775914d)

## DVGMRF
DVGMRFにおけるバイアスベクトルをガウス分布に従うと仮定した上で拡張したモデルです。
![dvhgmrf](https://github.com/tetdarth/image-denoising-for-GMRF/assets/136053901/aa7248e4-450a-4a72-b084-09b5286eba0a)


# ガウス型マルコフ確率場を用いた１次元信号ノイズ除去
## IVODGMRF
IVGMRFを1次元信号に適用したモデルです。
## DVODGMRF
DVGMRFを1次元信号に適用したモデルです。
## IVHODGMRF
IVHGMRFを1次元信号に適用したモデルです。
## DVHODGMRF
DVHGMRFを1次元信号に適用したモデルです。
