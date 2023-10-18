#include <fstream.h>
#include <string>
#include <vector>

int main()
{
    ifstream f_in("data/adult_tab.csv");
    ofstream f_out("data/adult_tab_numeric.txt");

    int age, fnlwgt, education_num, capital_gain, capital_loss, hours_per_week;
    string workclass, education, marital_status,
    occupation, relationship, race, sex, 
    native_country, income;

    while(f_in>>age>>workclass>>fnlwgt>>education>>education_num>>marital_status>>
        occupation>>relationship>>race>>sex>>capital_gain>>capital_loss>>hours_per_week
        native_country>>income)
    {
        f_out<<age<<"\t";
        int workclass_numeric = 0;

        if(workclass=="Private") 
            workclass_numeric = 1; 
        else
        if(workclass=="Self-emp-not-inc")
            workclass_numeric = 2; 
        else
        if(workclass=="Self-emp-inc") 
            workclass_numeric = 3; 
        else
        if(workclass=="Federal-gov")
            workclass_numeric = 4; 
        else
        if(workclass=="Local-gov")
            workclass_numeric = 5; 
        else
        if(workclass=="State-gov") 
            workclass_numeric = 6; 
        else
        if(workclass=="Without-pay")
            workclass_numeric = 7; 
        else
        if(workclass=="Never-worked")
            workclass_numeric = 8; 

        f_out<<workclass_numeric<<"\t"<<fnlwgt<<"\t";

        int education_numeric;

        if(education=="Bachelors")
            education_numeric=1;
        else
        if(education=="Some-college")
            education_numeric=2;
        else
        if(education=="11th")
            education_numeric=3;
        else
        if(education=="HS-grad")
            education_numeric=4;
        else
        if(education=="Prof-school")
            education_numeric=5;
        else
        if(education=="Assoc-acdm")
            education_numeric=6;
        else
        if(education=="Assoc-voc")
            education_numeric=7;
        else
        if(education=="9th")
            education_numeric=8;
        else
        if(education=="7th-8th")
            education_numeric=9;
        else
        if(education=="12th")
            education_numeric=10;
        else
        if(education=="Masters")
            education_numeric=11;
        else
        if(education=="1st-4th")
            education_numeric=12;
        else
        if(education=="10th")
            education_numeric=13;
        else
        if(education=="Doctorate")
            education_numeric=14;
        else
        if(education=="5th-6th")
            education_numeric=15;
        else
        if(education=="Preschool")
            education_numeric=16;

        f_out<<education_numeric<<"\t"<<education_num<<"\t";

        /*
marital-status	Feature	Categorical	Other	Married-civ-spouse, Divorced, Never-married, Separated, Widowed, Married-spouse-absent, Married-AF-spouse.		no
occupation	Feature	Categorical	Other	Tech-support, Craft-repair, Other-service, Sales, Exec-managerial, Prof-specialty, Handlers-cleaners, Machine-op-inspct, Adm-clerical, Farming-fishing, Transport-moving, Priv-house-serv, Protective-serv, Armed-Forces.		yes
relationship	Feature	Categorical	Other	Wife, Own-child, Husband, Not-in-family, Other-relative, Unmarried.		no
race	Feature	Categorical	Race	White, Asian-Pac-Islander, Amer-Indian-Eskimo, Other, Black.		no
sex	Feature	Binary	Sex	Female, Male.		no
capital-gain	Feature	Integer				no
capital-loss	Feature	Integer				no
hours-per-week	Feature	Integer				no
native-country	Feature	Categorical	Other	United-States, Cambodia, England, Puerto-Rico, Canada, Germany, Outlying-US(Guam-USVI-etc), India, Japan, Greece, South, China, Cuba, Iran, Honduras, Philippines, Italy, Poland, Jamaica, Vietnam, Mexico, Portugal, Ireland, France, Dominican-Republic, Laos, Ecuador, Taiwan, Haiti, Columbia, Hungary, Guatemala, Nicaragua, Scotland, Thailand, Yugoslavia, El-Salvador, Trinadad&Tobago, Peru, Hong, Holand-Netherlands.		yes
income	Target	Binary	Income	>50K, <=50K.		no        
        */

        f_out<<endl;
    }

}