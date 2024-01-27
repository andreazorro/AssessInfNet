# AssessInfNet
Assessment of inferred networks in Matlab

Required Input 
--------------

* One or multiple *.tsv files with the inferred interactions. Only first two columns are taking into account. File name "Method_data.tsv" (See file Anova_Scoelicolor-rmabatch.tsv)
* One or multiple *.gs files with the gold standard for the assessment. Only first two columns are taking into account. 

Output
------

* A *.txt file with the results of the assesment (TP: True Positives, Recall, MCC and F1 metrics)
* One or multiple heatmaps of the metrics for each gold standard. 

License
-------

This project is licensed under the GNU General Public License. For the exact terms please see the [LICENSE file]
