import pandas as pd
from sdv.datasets.local import load_csvs
from sdv.metadata import SingleTableMetadata
from sdv.sequential import PARSynthesizer

datasets = load_csvs(folder_name="../data/")
data = datasets["liver_cirrhosis"]
data_pd = pd.read_csv("../data/liver_cirrhosis.csv")
metadata = SingleTableMetadata()

metadata.detect_from_dataframe(data=data_pd)
metadata.update_column(column_name="patient_id", sdtype="id")

metadata.set_sequence_key(column_name="patient_id")
metadata.set_sequence_index(column_name="Tstop")

synthesizer = PARSynthesizer(
    metadata, context_columns=["age", "female", "treatment_grp", "starting_state"]
)
synthesizer.fit(data_pd)
synthesizer.save(filepath="../models/saved_synthesizer_sdvseq.pkl")

synthetic_data = synthesizer.sample(num_sequences=100)
synthetic_data.to_csv("../data/synthetic_data_sdvseq.csv", index=False)
