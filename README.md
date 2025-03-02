# Differential diagnosis based on a patient's symptoms

## Project Overview

Differential diagnosis based on a patient's symptoms is a Streamlit-based web application that facilitates phenotype-driven disease matching using the Human Phenotype Ontology (HPO). This tool allows healthcare professionals and researchers to input patient phenotypes and match them against known genetic disorders, providing valuable insights for diagnosis and research.

## Features

- Clinical notes to HPO term conversion
- Manual HPO term entry
- Automated disease matching based on patient phenotypes
- AI-powered symptom recommendations
- Visualization of disease networks and HPO structure
- Real-time similarity scoring

## Technologies Used

- Python 3.8+
- Streamlit
- Requests
- JSON
- Regular Expressions (re)
- Custom PhrankCalculator class
- Local Llama API integration

## External APIs and Services

- doc2hpo API for HPO term extraction
- Local Llama API for AI-powered recommendations

## Data Files

- HPO OBO file (hp.obo)
- HPO Annotation file (phenotype.hpoa)

## Setup and Installation

1. Clone the repository
2. Install required packages: `pip install streamlit requests`
3. Ensure you have the necessary data files (hp.obo and phenotype.hpoa) in the `data/` directory
4. Set up and run the local Llama API server

## Running the Application

Execute the following command in the project directory:

```
streamlit run app.py
```

## Main Components

1. **Patient Input Section**: Allows users to input clinical notes or manually enter HPO terms.
2. **HPO Term Generation**: Converts clinical notes to HPO terms using external APIs.
3. **Disease Matching**: Utilizes the PhrankCalculator to rank diseases based on patient phenotypes.
4. **AI Recommendations**: Integrates with Llama API to suggest additional relevant symptoms.
5. **Visualization**: Provides network and structure visualizations of diseases and HPO terms.
6. **Analysis Results**: Displays top disease matches and similarity scores.

## Note

This application requires specific data files and API access. Ensure all dependencies are correctly set up before running the application.

Sources
[1] Build a Python Website in 15 Minutes With Streamlit - YouTube https://www.youtube.com/watch?v=2siBrMsqF44
[2] Python Requests Module - W3Schools https://www.w3schools.com/python/module_requests.asp
[3] Python Json Module - Python Cheatsheet https://www.pythoncheatsheet.org/modules/json-module
[4] Python Regular Expressions re.match() and re.sub() Explained https://builtin.com/articles/python-re-match
[5] Class Rank Calculator | Free Academic Standing & Percentile Tool https://classrankcalculator.xyz
[6] About Streamlit in Snowflake https://docs.snowflake.com/en/developer-guide/streamlit/about-streamlit
[7] Requests - PyPI https://pypi.org/project/requests/
[8] Python JSON Data: A Guide With Examples - DataCamp https://www.datacamp.com/tutorial/json-data-python
[9] re – Regular Expressions - Python Module of the Week - PyMOTW 3 https://pymotw.com/2/re/
[10] Density Calculator - Determine Your Freight Class https://shipgsl.com/shipping-tools/density-calculator/
[11] Streamlit • A faster way to build and share data apps https://streamlit.io
[12] Requests: HTTP for Humans™ — Requests 2.32.3 documentation https://requests.readthedocs.io
[13] Python JSON - W3Schools https://www.w3schools.com/python/python_json.asp
[14] Python Regular Expressions - Google for Developers https://developers.google.com/edu/python/regular-expressions
[15] Finding Probabilities for Normally Distributed Data Using ClassCalc https://www.youtube.com/watch?v=WenUs0sCuho
[16] Streamlit — A faster way to build and share data apps. - GitHub https://github.com/streamlit/streamlit
[17] I don't really understand Requests python module - Reddit https://www.reddit.com/r/learnpython/comments/w6q566/i_dont_really_understand_requests_python_module/
[18] json — JSON encoder and decoder — Python 3.13.2 documentation https://docs.python.org/3/library/json.html
[19] re — Regular expression operations — Python 3.13.2 documentation https://docs.python.org/3/library/re.html
[20] ClassCalc - Test Safe Online Graphing Calculator https://classcalc.com
