python CPPred.py -i ../data/Human_coding_RNA_test.fa  -hex ../Hexamer/Human_Hexamer.tsv -r ../Human_Model/Human.range -mol ../Human_Model/Human.model -spe Human -o Human_coding.result 
python CPPred.py -i ../data/Human_ncrna_test.fa -hex ../Hexamer/Human_Hexamer.tsv -r ../Human_Model/Human.range -mol ../Human_Model/Human.model -spe Human -o Human_ncrna.result 

python CPPred.py -i ../data/Mouse_coding_RNA_test.fa  -hex ../Hexamer/Human_Hexamer.tsv -r ../Human_Model/Human.range -mol ../Human_Model/Human.model -spe Human -o Mouse_coding.result 
python CPPred.py -i ../data/Mouse_ncrna_test.fa -hex ../Hexamer/Human_Hexamer.tsv -r ../Human_Model/Human.range -mol ../Human_Model/Human.model -spe Human -o Mouse_ncrna.result 

python CPPred.py -i ../data/Zebrasfish_coding_RNA_test.fa  -hex ../Hexamer/Human_Hexamer.tsv -r ../Human_Model/Human.range -mol ../Human_Model/Human.model -spe Human -o Zebrafish_coding.result 
python CPPred.py -i ../data/Zebrafish_ncrna_test.fa -hex ../Hexamer/Human_Hexamer.tsv -r ../Human_Model/Human.range -mol ../Human_Model/Human.model -spe Human -o Zebrafish_ncrna.result 

python CPPred.py -i ../data/S.cerevisiae_coding_RNA_test.fa  -hex ../Hexamer/Human_Hexamer.tsv -r ../Human_Model/Human.range -mol ../Human_Model/Human.model -spe Human -o S.cerevisiae_coding.result 
python CPPred.py -i ../data/S.cerevisiae_ncrna_test.fa -hex ../Hexamer/Human_Hexamer.tsv -r ../Human_Model/Human.range -mol ../Human_Model/Human.model -spe Human -o S.cerevisiae_ncrna.result 

python CPPred.py -i ../data/Fruit_fly_coding_RNA_test.fa  -hex ../Hexamer/Human_Hexamer.tsv -r ../Human_Model/Human.range -mol ../Human_Model/Human.model -spe Human -o Fruit_fly_coding.result 
python CPPred.py -i ../data/Fruit_fly_ncrna_test.fa -hex ../Hexamer/Human_Hexamer.tsv -r ../Human_Model/Human.range -mol ../Human_Model/Human.model -spe Human -o Fruit_fly_ncrna.result 

python CPPred.py -i ../data/Integrated_coding_RNA_test.fa  -hex ../Hexamer/Integrated_Hexamer.tsv -r ../Integrated_Model/Integrated.range -mol ../Integrated_Model/Integrated.model -spe Integrated -o Integrated_coding.result 
python CPPred.py -i ../data/Integrated_ncrna_test.fa -hex ../Hexamer/Integrated_Hexamer.tsv -r ../Integrated_Model/Integrated.range -mol ../Integrated_Model/Integrated.model -spe Integrated -o Integrated_ncrna.result 
