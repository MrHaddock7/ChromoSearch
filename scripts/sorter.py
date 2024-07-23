import pandas as pd
import logging

logger = logging.getLogger(__name__)

# Load the CSV file into a DataFrame
# input_file = 'blastp_results.csv'  # Replace with your input file path
# output_file = 'sorted_output_names.csv'  # Replace with your desired output file path

def csv_sorter(input_csv, genome, output, shortest_sequence=0, e_value=0.05):
    logger.debug(f'Entering csv_sorter function')
    try:
        # Read the CSV file
        df = pd.read_csv(input_csv)

        filtered_df = df[df['evalue'] < e_value]
        df_sorted = filtered_df.sort_values(by='evalue', ascending=False)

        # Write the sorted DataFrame to a new CSV file
        df_sorted.to_csv(f'{output}/output_{genome}_sorter.csv', index=False)

        print(f"CSV file sorted by column 'e' and saved to output_{genome}_sorter.csv")
        logger.info(f'Results from csv_sorter function are saved in file: output_{genome}_sorter.csv')
    except Exception as e:
        logger.error(f'Error in csv_sorter function: {e}')
    except KeyboardInterrupt:
        logger.warning("Data processing interrupted by user")
    logger.debug(f'Exiting csv_sorter function')

def csv_sorter2(input_csv, genome, output):
    logger.debug(f'Entering csv_sorter2 function')
    try:
        # Read the CSV file
        df = pd.read_csv(input_csv)

        df_sorted = df.sort_values(by='Score', ascending=False)

        # Write the sorted DataFrame to a new CSV file
        df_sorted.to_csv(f'{output}/output_{genome}_sorter2.csv', index=False)

        print(f"CSV file sorted by column 'e' and saved to output_{genome}_sorter.csv")
        logger.info(f'Results from csv_sorter function are saved in file: output_{genome}_sorter2.csv')
    except Exception as ex:
        logger.error(f'Error in csv_sorter2 function: {ex}')
    except KeyboardInterrupt:
        logger.warning("Data processing interrupted by user")
    logger.debug(f'Exiting csv_sorter2 function')