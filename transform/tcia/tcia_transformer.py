from transform.tcia.tcia_harvester import load_all
import networkx as nx


class TCIATransformer():

    def __init__(self):
        # fetch data
        (self.all_patients, self.all_patient_study, self.all_series) = load_all()
        for a in (self.all_patients, self.all_patient_study, self.all_series):
            print(len(a))

    def graph(self):
        """Creates graph"""
        G = nx.MultiDiGraph()
        G.add_node('tcia', label= 'Source')
        
        G.add_nodes_from(set([case.Collection for case in self.all_patients]), label='Project')
        G.add_edges_from([('tcia',k) for k in set([case.Collection for case in self.all_patients])], label='projects')

        G.add_nodes_from(set([case.PatientID for case in self.all_patients]), label='Case')
        G.add_edges_from([(case.Collection, case.PatientID) for case in self.all_patients], label='cases')

        # add ['PatientName', 'PatientSex']
        G.add_nodes_from(set([case.PatientID + '-demographic' for case in self.all_patients]), label='Demographic')
        G.add_edges_from([(case.PatientID, case.PatientID + '-demographic') for case in self.all_patients], label='demographics')

        G.add_nodes_from(set([study.StudyInstanceUID for study in self.all_patient_study]), label='Study')
        G.add_edges_from([(study.PatientID, study.StudyInstanceUID) for study in self.all_patient_study], label='studies')

        # note the api seems to be returning one series that is really a study
        series = list(filter(lambda s: s.SeriesInstanceUID != '1.3.6.1.4.1.14519.5.2.1.6279.6001.287560874054243719452635194040', self.all_series))
        G.add_nodes_from(set([s.SeriesInstanceUID for s in series]), label='Series')
        G.add_edges_from([(s.StudyInstanceUID, s.SeriesInstanceUID) for s in series], label='series')

        return G
