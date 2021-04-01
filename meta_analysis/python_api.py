#/Users/dutty/.pyenv/shims/python

from serpapi import GoogleSearch

params = {
  "api_key": "8ead96006c7d5488de24656bf2224b843c86208cdb8589258b39cf9bfacee1db",
  "engine": "google_scholar",
  "as_ylo": "2000",
  "as_yhi": "2021",
  "q": ""balancing selection"",
  "hl": "en"
}

search = GoogleSearch(params)
results = search.get_dict()
