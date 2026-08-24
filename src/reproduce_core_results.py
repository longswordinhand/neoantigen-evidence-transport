#!/usr/bin/env python3
"""Reproduce the core cross-trial transportability analyses from the harmonized table.

The script intentionally starts from the public derived benchmark rather than redistributing
third-party supplementary files. It performs leave-one-trial-out prediction, within-patient
ranking, trial-identifiability analysis, peptide-length ablation, and the study-native
epitope-coverage sensitivity analysis used in the manuscript.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.ensemble import ExtraTreesClassifier, RandomForestClassifier, HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, brier_score_loss, log_loss, roc_auc_score
from sklearn.model_selection import GroupKFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "harmonized_vaccine_targets.csv"
OUT = ROOT / "results" / "reproduced_summary.json"
SEED = 17

D = pd.read_csv(DATA)
TRIALS = list(D["trial"].unique())
Y = D["response"].astype(int).to_numpy()
COMPOSITION = [c for c in D.columns if c.startswith("frac_")] + ["aa_entropy"]
FULL_SEQUENCE = ["length"] + COMPOSITION


def lr_model():
    return Pipeline([
        ("scale", StandardScaler()),
        ("model", LogisticRegression(C=1.0, max_iter=5000, class_weight="balanced", solver="liblinear", random_state=SEED)),
    ])


def metrics(y, p):
    return {
        "n": int(len(y)),
        "positive": int(np.sum(y)),
        "auc": float(roc_auc_score(y, p)),
        "logloss": float(log_loss(y, p, labels=[0, 1])),
        "brier": float(brier_score_loss(y, p)),
    }


def loto_predictions(features, estimator):
    pred = np.full(len(D), np.nan)
    per_trial = {}
    for held in TRIALS:
        te = D["trial"].eq(held).to_numpy()
        tr = ~te
        Xtr = D.loc[tr, features].copy()
        Xte = D.loc[te, features].copy()
        med = Xtr.median(numeric_only=True)
        Xtr = Xtr.fillna(med).fillna(0)
        Xte = Xte.fillna(med).fillna(0)
        model = clone(estimator)
        model.fit(Xtr, Y[tr])
        p = model.predict_proba(Xte)[:, 1]
        pred[te] = p
        per_trial[held] = metrics(Y[te], p)
    return pred, per_trial


def macro_patient_auc(pred):
    result = {}
    for trial in TRIALS:
        z = D.loc[D.trial.eq(trial), ["patient", "response"]].copy()
        z["pred"] = pred[D.trial.eq(trial).to_numpy()]
        aucs = []
        for _, g in z.groupby("patient"):
            if g.response.nunique() < 2:
                continue
            aucs.append(roc_auc_score(g.response, g.pred))
        result[trial] = {"eligible_patients": len(aucs), "macro_patient_auc": float(np.mean(aucs)) if aucs else None}
    return result


def trial_identity(features):
    X = D[features].fillna(0.0)
    y = D.trial.astype(str)
    groups = D.trial.astype(str) + "::" + D.patient.astype(str)
    model = Pipeline([("scale", StandardScaler()), ("model", LogisticRegression(max_iter=5000, random_state=SEED))])
    p = cross_val_predict(model, X, y, groups=groups, cv=GroupKFold(n_splits=5), method="predict")
    return float(balanced_accuracy_score(y, p))


out = {"seed": SEED, "n_targets": int(len(D)), "n_patients": int((D.trial.astype(str) + "::" + D.patient.astype(str)).nunique())}

for label, features in {"full_sequence": FULL_SEQUENCE, "composition_no_length": COMPOSITION}.items():
    pred, fold = loto_predictions(features, lr_model())
    out[label] = {
        "leave_one_trial_out": fold,
        "macro_within_patient": macro_patient_auc(pred),
        "pooled_auc": float(roc_auc_score(Y, pred)),
    }

out["trial_identity_balanced_accuracy"] = {
    "length_only": trial_identity(["length"]),
    "composition_no_length": trial_identity(COMPOSITION),
    "full_sequence": trial_identity(FULL_SEQUENCE),
    "chance_three_class": 1.0 / 3.0,
}

models = {
    "logistic": lr_model(),
    "random_forest": RandomForestClassifier(n_estimators=500, min_samples_leaf=3, class_weight="balanced", random_state=SEED, n_jobs=-1),
    "extra_trees": ExtraTreesClassifier(n_estimators=500, min_samples_leaf=3, class_weight="balanced", random_state=SEED, n_jobs=-1),
    "hist_gradient_boosting": HistGradientBoostingClassifier(random_state=SEED),
}
out["model_stability"] = {}
for name, model in models.items():
    pred, fold = loto_predictions(COMPOSITION, model)
    mp = macro_patient_auc(pred)
    out["model_stability"][name] = {
        "mean_trial_auc": float(np.mean([v["auc"] for v in fold.values()])),
        "mean_macro_patient_auc": float(np.mean([v["macro_patient_auc"] for v in mp.values()])),
        "pooled_auc": float(roc_auc_score(Y, pred)),
        "per_trial": fold,
    }

for label, features in {"composition": COMPOSITION, "composition_plus_epitope_coverage": COMPOSITION + ["epitope_coverage"]}.items():
    pred, fold = loto_predictions(features, lr_model())
    out.setdefault("epitope_coverage_sensitivity", {})[label] = {
        "leave_one_trial_out": fold,
        "macro_within_patient": macro_patient_auc(pred),
    }

OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(json.dumps(out, indent=2) + "\n")
print(json.dumps(out, indent=2))
