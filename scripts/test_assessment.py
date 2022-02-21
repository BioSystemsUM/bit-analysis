from unittest import TestCase

from scripts.models_comparison import getMtuberculosis_models
from sklearn.metrics import balanced_accuracy_score, precision_score, recall_score, matthews_corrcoef
from assessment import ResultsReport


class TestAssessment(TestCase):

    def test_radar_graph(self):
        Mtuberculosis_models = getMtuberculosis_models()

        metrics = [balanced_accuracy_score, precision_score, recall_score, matthews_corrcoef]
        report = ResultsReport(Mtuberculosis_models["reference model"],
                               [Mtuberculosis_models["carveme model"],
                                Mtuberculosis_models["all model"]])

        report.generate_radar_graph(metrics, ["Balanced Accuracy", "Precision", "Recall", "MCC"], "radar.png")
