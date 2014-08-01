echo 'actPro'
python computePostPr.py AtacPkTag_actPro.bed AtacPkTag_hetero.bed AtacPkTag_actPro.bed |cut -f4 |awk 'BEGIN{printf "actPro"}{printf "\t"$1}END{printf "\n"}' >temp.dat

echo 'hetero'
python computePostPr.py AtacPkTag_actPro.bed AtacPkTag_hetero.bed AtacPkTag_hetero.bed |cut -f4 |awk 'BEGIN{printf "hetero"}{printf "\t"$1}END{printf "\n"}' >>temp.dat

R --no-save --slave --args temp.dat <~/bin/plotBox.R
