

shots = [79797, 80243, 80257, 80328, 80376, 80745]

from DBS.io.DBSsync import push_processed_data
for shot in shots:
    push_processed_data('tcv', shot, data='beamtracing', dest_hostname='altair1')
    push_processed_data('tcv', shot, data='fDop_estimation', dest_hostname='altair1')
    
    