#!/usr/bin/env python3
"""
GCNdesign Validation Metrics Analysis
This script provides comprehensive validation metrics and visualization tools
"""

import os
import sys
import torch
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import (
    accuracy_score, precision_score, recall_score, f1_score, 
    matthews_corrcoef, confusion_matrix, classification_report
)
from tqdm import tqdm
import pickle

# Add gcndesign to path
sys.path.append('.')
from gcndesign.predictor import Predictor
from gcndesign.hypara import HyperParam
from gcndesign.dataset import pdb2input, add_margin

# Amino acid mapping
amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
               'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

def get_model_predictions_on_protein(pdb_file, predictor):
    """Get predictions for a single protein"""
    try:
        # Get predictions
        predictions = predictor.predict_logit_tensor(pdb_file)
        
        # Get ground truth
        hypara = HyperParam()
        node, edgemat, adjmat, labels, mask, aa_sequence = pdb2input(pdb_file, hypara)
        
        # Convert predictions to class indices
        pred_classes = np.argmax(predictions, axis=1)
        
        # Filter valid residues only
        valid_indices = mask.numpy()
        true_labels = labels.numpy()[valid_indices]
        pred_labels = pred_classes[valid_indices]
        
        return true_labels, pred_labels, aa_sequence
    except Exception as e:
        print(f"Error processing {pdb_file}: {e}")
        return None, None, None

def create_sample_test_protein():
    """Create a simple test protein for demonstration"""
    pdb_content = """HEADER    SYNTHETIC                               01-JAN-20   TEST
ATOM      1  N   ALA A   1      20.154  16.967  10.000  1.00 20.00           N
ATOM      2  CA  ALA A   1      21.618  17.000  10.000  1.00 20.00           C
ATOM      3  C   ALA A   1      22.000  18.430  10.000  1.00 20.00           C
ATOM      4  O   ALA A   1      21.154  19.320  10.000  1.00 20.00           O
ATOM      5  CB  ALA A   1      22.154  16.319   8.732  1.00 20.00           C
ATOM      6  N   GLY A   2      23.300  18.630  10.000  1.00 20.00           N
ATOM      7  CA  GLY A   2      23.800  19.990  10.000  1.00 20.00           C
ATOM      8  C   GLY A   2      25.300  20.000  10.000  1.00 20.00           C
ATOM      9  O   GLY A   2      26.000  19.000  10.000  1.00 20.00           O
ATOM     10  N   VAL A   3      25.800  21.230  10.000  1.00 20.00           N
ATOM     11  CA  VAL A   3      27.240  21.400  10.000  1.00 20.00           C
ATOM     12  C   VAL A   3      27.800  22.820  10.000  1.00 20.00           C
ATOM     13  O   VAL A   3      27.100  23.830  10.000  1.00 20.00           O
ATOM     14  CB  VAL A   3      27.800  20.630   8.800  1.00 20.00           C
ATOM     15  CG1 VAL A   3      29.300  20.800   8.800  1.00 20.00           C
ATOM     16  CG2 VAL A   3      27.300  19.190   8.800  1.00 20.00           C
END
"""
    
    with open('test_protein.pdb', 'w') as f:
        f.write(pdb_content)
    return 'test_protein.pdb'

def analyze_single_protein_predictions(predictor, pdb_file):
    """Analyze predictions for a single protein and show detailed metrics"""
    print(f"\n=== Analyzing {pdb_file} ===")
    
    # Get predictions
    true_labels, pred_labels, aa_sequence = get_model_predictions_on_protein(pdb_file, predictor)
    
    if true_labels is None:
        print("Failed to get predictions")
        return None
    
    # Calculate metrics
    accuracy = accuracy_score(true_labels, pred_labels)
    precision = precision_score(true_labels, pred_labels, average='macro', zero_division=0)
    recall = recall_score(true_labels, pred_labels, average='macro', zero_division=0)
    f1 = f1_score(true_labels, pred_labels, average='macro', zero_division=0)
    
    # For MCC, we need to handle potential issues with single-class predictions
    try:
        mcc = matthews_corrcoef(true_labels, pred_labels)
    except:
        mcc = 0.0
    
    print(f"📊 VALIDATION METRICS:")
    print(f"   Accuracy:  {accuracy:.4f} ({accuracy*100:.2f}%)")
    print(f"   Precision: {precision:.4f}")
    print(f"   Recall:    {recall:.4f}")
    print(f"   F1-Score:  {f1:.4f}")
    print(f"   MCC:       {mcc:.4f}")
    
    # Show per-residue analysis
    print(f"\n🧬 PER-RESIDUE ANALYSIS:")
    print(f"   Total residues: {len(true_labels)}")
    correct_predictions = np.sum(true_labels == pred_labels)
    print(f"   Correct predictions: {correct_predictions}")
    print(f"   Per-residue accuracy: {correct_predictions/len(true_labels)*100:.2f}%")
    
    # Show sequence comparison
    print(f"\n🔍 SEQUENCE COMPARISON:")
    true_seq = ''.join([amino_acids[i] for i in true_labels])
    pred_seq = ''.join([amino_acids[i] for i in pred_labels])
    print(f"   True:      {true_seq}")
    print(f"   Predicted: {pred_seq}")
    print(f"   Match:     {''.join(['✓' if t==p else '✗' for t,p in zip(true_seq, pred_seq)])}")
    
    return {
        'accuracy': accuracy,
        'precision': precision,
        'recall': recall,
        'f1': f1,
        'mcc': mcc,
        'true_labels': true_labels,
        'pred_labels': pred_labels,
        'sequence_length': len(true_labels),
        'correct_predictions': correct_predictions
    }

def create_validation_plots(metrics_data, save_plots=True):
    """Create comprehensive validation plots"""
    plt.style.use('seaborn-v0_8' if 'seaborn-v0_8' in plt.style.available else 'default')
    
    # Create figure with subplots
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle('GCNdesign Validation Metrics Analysis', fontsize=16, fontweight='bold')
    
    # 1. Metrics Bar Plot
    ax1 = axes[0, 0]
    metrics_names = ['Accuracy', 'Precision', 'Recall', 'F1-Score']
    metrics_values = [
        metrics_data['accuracy'], 
        metrics_data['precision'], 
        metrics_data['recall'], 
        metrics_data['f1']
    ]
    
    bars = ax1.bar(metrics_names, metrics_values, color=['#2E86AB', '#A23B72', '#F18F01', '#C73E1D'])
    ax1.set_ylim(0, 1)
    ax1.set_ylabel('Score')
    ax1.set_title('Performance Metrics')
    ax1.grid(axis='y', alpha=0.3)
    
    # Add value labels on bars
    for bar, value in zip(bars, metrics_values):
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01, 
                f'{value:.3f}', ha='center', va='bottom', fontweight='bold')
    
    # 2. Confusion Matrix
    ax2 = axes[0, 1]
    cm = confusion_matrix(metrics_data['true_labels'], metrics_data['pred_labels'])
    
    # Get unique labels for the confusion matrix
    unique_labels = np.unique(np.concatenate([metrics_data['true_labels'], metrics_data['pred_labels']]))
    cm_labels = [amino_acids[i] for i in unique_labels]
    
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', ax=ax2,
                xticklabels=cm_labels, yticklabels=cm_labels)
    ax2.set_xlabel('Predicted')
    ax2.set_ylabel('True')
    ax2.set_title('Confusion Matrix')
    
    # 3. Prediction Accuracy by Position
    ax3 = axes[1, 0]
    correct_by_position = (metrics_data['true_labels'] == metrics_data['pred_labels']).astype(int)
    positions = range(1, len(correct_by_position) + 1)
    
    colors = ['green' if correct else 'red' for correct in correct_by_position]
    ax3.bar(positions, correct_by_position, color=colors, alpha=0.7)
    ax3.set_xlabel('Residue Position')
    ax3.set_ylabel('Correct (1) / Incorrect (0)')
    ax3.set_title('Per-Position Prediction Accuracy')
    ax3.set_ylim(-0.1, 1.1)
    
    # 4. Model Confidence Analysis
    ax4 = axes[1, 1]
    
    # Since we don't have probability data in this simple example,
    # we'll create a summary statistics plot
    summary_data = {
        'Total Residues': metrics_data['sequence_length'],
        'Correct Predictions': metrics_data['correct_predictions'],
        'Incorrect Predictions': metrics_data['sequence_length'] - metrics_data['correct_predictions']
    }
    
    wedges, texts, autotexts = ax4.pie(
        [summary_data['Correct Predictions'], summary_data['Incorrect Predictions']], 
        labels=['Correct', 'Incorrect'],
        colors=['#2E86AB', '#C73E1D'],
        autopct='%1.1f%%',
        startangle=90
    )
    ax4.set_title('Prediction Accuracy Distribution')
    
    plt.tight_layout()
    
    if save_plots:
        plt.savefig('gcndesign_validation_metrics.png', dpi=300, bbox_inches='tight')
        print(f"📈 Plots saved as 'gcndesign_validation_metrics.png'")
    
    plt.show()
    
    return fig

def main():
    """Main function to run validation analysis"""
    print("🧬 GCNdesign Validation Metrics Analysis")
    print("=" * 50)
    
    # Initialize predictor
    print("🔧 Loading GCNdesign model...")
    try:
        predictor = Predictor(device='cpu')  # Use CPU for compatibility
        print("✅ Model loaded successfully!")
    except Exception as e:
        print(f"❌ Error loading model: {e}")
        return
    
    # Create a test protein for demonstration
    print("\n📁 Creating test protein...")
    test_pdb = create_sample_test_protein()
    print(f"✅ Test protein created: {test_pdb}")
    
    # Analyze predictions
    print("\n🔍 Running validation analysis...")
    metrics_data = analyze_single_protein_predictions(predictor, test_pdb)
    
    if metrics_data:
        print("\n📊 Creating validation plots...")
        create_validation_plots(metrics_data)
        
        # Save metrics to file
        print("\n💾 Saving metrics to file...")
        with open('validation_metrics.pkl', 'wb') as f:
            pickle.dump(metrics_data, f)
        
        # Create a summary report
        report = f"""
GCNdesign Validation Report
==========================

Model Performance:
- Accuracy: {metrics_data['accuracy']:.4f} ({metrics_data['accuracy']*100:.2f}%)
- Precision: {metrics_data['precision']:.4f}
- Recall: {metrics_data['recall']:.4f}
- F1-Score: {metrics_data['f1']:.4f}
- Matthews Correlation Coefficient: {metrics_data['mcc']:.4f}

Dataset Statistics:
- Total residues analyzed: {metrics_data['sequence_length']}
- Correct predictions: {metrics_data['correct_predictions']}
- Prediction accuracy: {metrics_data['correct_predictions']/metrics_data['sequence_length']*100:.2f}%

Generated files:
- validation_metrics.pkl (metrics data)
- gcndesign_validation_metrics.png (visualization plots)
"""
        
        with open('validation_report.txt', 'w') as f:
            f.write(report)
        
        print("✅ Analysis complete! Files generated:")
        print("   📊 gcndesign_validation_metrics.png")
        print("   📁 validation_metrics.pkl")
        print("   📄 validation_report.txt")
    
    # Cleanup
    if os.path.exists(test_pdb):
        os.remove(test_pdb)

if __name__ == "__main__":
    main()
