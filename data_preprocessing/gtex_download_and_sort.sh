## run this in the directory you want your eQTLs to be in

wget http://www.gtexportal.org/static/datasets/gtex_analysis_v6/single_tissue_eqtl_data/GTEx_Analysis_V6_eQTLs.tar.gz
tar -xzvf GTEx_Analysis_V6_eQTLs.tar.gz
mkdir sorted
for f in *.snpgenes; do
    FNAME=`basename $f`
    echo "Sorting file ${FNAME}"
    sort -k14,14n -k15,15n $f > sorted/$FNAME
    rm $f
done
