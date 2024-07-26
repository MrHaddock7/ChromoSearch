import pandas as pd

# Sample DataFrames
df1 = pd.DataFrame({
    'Protein Name': [
        'sp|P00914|PHR_ECOLI', 
        'sp|P05066|PHR_YEAST', 
        'sp|P06585|PSBA_PEA', 
        'sp|P0A6D3|AROA_ECOLI', 
        'sp|P0A7W1|RS5_ECOLI'
    ],
    'Sequence': [
        'MTTHLVWFRQDLRLHDNLALAAACRNSSARVLALYIATPRQWATHN...',
        'MKRTVISSSNAYASKRSRLDIEHDFEQYHSLNKKYYPRPITRTGAN...',
        'MTAILERRDSENLWGRFCNWITSTENRLYIGWFGVLMIPTLLTATS...',
        'MESLTLQPIARVDGTINLPGSKSVSNRALLLAALAHGKTVLTNLLD...',
        'MAHIEKQAGELQEKLIAVNRVSKTVKGGRIFSFTALTVVGDGNGRV...'
    ]
})

df2 = pd.DataFrame({
    'Protein Name': [
        'sp|P12345|PHR_HUMAN', 
        'sp|P23456|PHR_MOUSE', 
        'sp|P34567|PHR_RAT', 
        'sp|P45678|PHR_FISH', 
        'sp|P56789|PHR_BIRD'
    ],
    'Sequence': [
        'MKKLLKWFRLDSLDLDNLALAAACRNSSARVLALYIATPRQWATHN...',
        'MKRTLISSNAYASKRSRLDIEHDFEQYHSLNKKYYPRPITRTGAN...',
        'MTAILERDSENLWGRFCNWITSTENRLYIGWFGVLMIPTLLTATS...',
        'MESLTQPIARVDGTINLPGSKSVSNRALLLAALAHGKTVLTNLLD...',
        'MAHIEKQAGELQEKLIAVNRVSKTVKGGRIFSFTALTVVGDGNGRV...'
    ]
})

# Using pd.merge with how='cross' to get the Cartesian product
cross_product_df = pd.merge(df1, df2, how='cross')

# Renaming columns for clarity
cross_product_df.columns = [
    'Protein Name 1', 'Sequence 1',
    'Protein Name 2', 'Sequence 2'
]

print(cross_product_df)
