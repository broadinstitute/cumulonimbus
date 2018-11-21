import googleapiclient.discovery

project_id = 'broad-gdr-dig'
region = 'Multi-Region'
zone = '???'

dataproc = googleapiclient.discovery.build('dataproc', 'v1')
