from pandas import DataFrame,read_csv

# http://data.humdata.org/dataset/798bb02a-3363-4267-9e4e-b03e7bd2195d/resource/8059ec97-41ba-4ea6-a300-d064fa91d8df/download/covid-19-historical-cases-by-country.csv


if __name__ == "__main__":

    df = read_csv(r'./covid-19-historical-cases-by-country.csv')
    df['ADM0_NAME']=df['ADM0_NAME'].str.upper()
    df[['ADM0_NAME','DateOfDataEntry','cum_conf']].rename(columns={'ADM0_NAME':'country','DateOfDataEntry':'date','cum_conf':'cases'}).sort_values(by=['country','date']).\
        to_csv(r'./humdata_case_numbers.csv',index=False,sep=';')
