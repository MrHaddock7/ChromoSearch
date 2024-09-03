import pandas as pd
import logging

logger = logging.getLogger(__name__)

def csv_sorter(input_csv, genome, output, sort_value_metric, name_output, cut_off_value=False, greater_than=False, only_sort=False):
    try:
        df = pd.read_csv(input_csv)
        df_sorted = df.sort_values(by=sort_value_metric, ascending=False)
        if not only_sort:
            if greater_than:
                df_sorted = df_sorted[df_sorted[sort_value_metric] > cut_off_value]
                df_sorted = df_sorted.sort_values(by=sort_value_metric, ascending=True)
            else:
                df_sorted = df_sorted[df_sorted[sort_value_metric] < cut_off_value]
                df_sorted = df_sorted.sort_values(by=sort_value_metric, ascending=True)
        logger.info(f'Results from csv_sorter function are saved in file: output_{genome}_{name_output}.csv')
    except Exception as ex:
        logger.error(f'Error in csv_sorter function: {ex}')
    df_sorted.to_csv(f'{output}/output_{genome}_{name_output}.csv', index=False)


      
    