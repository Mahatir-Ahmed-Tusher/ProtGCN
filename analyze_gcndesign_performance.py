#!/usr/bin/env python3
"""
GCNdesign Performance Analysis with Real Validation Metrics
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
import pickle

# Add gcndesign to path
sys.path.append('.')
from gcndesign.predictor import Predictor
from gcndesign.hypara import HyperParam
from gcndesign.dataset import pdb2input, add_margin

# Amino acid mapping
amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
               'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

def analyze_protein_predictions(pdb_file):
    """Comprehensive analysis of GCNdesign predictions"""
    print(f"\n🧬 ANALYZING: {pdb_file}")
    print("=" * 60)
    
    # Initialize predictor
    predictor = Predictor(device='cpu')
    hypara = HyperParam()
    
    # Get predictions and ground truth
    print("🔍 Getting model predictions...")
    predictions = predictor.predict_logit_tensor(pdb_file)
    pred_probs = torch.softmax(torch.tensor(predictions), dim=1).numpy()
    
    print("🔍 Getting ground truth labels...")
    node, edgemat, adjmat, labels, mask, aa_sequence = pdb2input(pdb_file, hypara)
    
    # Get valid predictions only
    valid_mask = mask.numpy()
    true_labels = labels.numpy()[valid_mask]
    pred_labels = np.argmax(predictions, axis=1)[valid_mask]
    pred_probs_valid = pred_probs[valid_mask]
    
    print(f"📊 Protein Statistics:")
    print(f"   Total residues: {len(aa_sequence)}")
    print(f"   Valid residues: {len(true_labels)}")
    print(f"   Sequence: {''.join(aa_sequence[:50])}{'...' if len(aa_sequence) > 50 else ''}")
    
    # Calculate metrics
    accuracy = accuracy_score(true_labels, pred_labels)
    precision = precision_score(true_labels, pred_labels, average='macro', zero_division=0)
    recall = recall_score(true_labels, pred_labels, average='macro', zero_division=0)
    f1 = f1_score(true_labels, pred_labels, average='macro', zero_division=0)
    
    try:
        mcc = matthews_corrcoef(true_labels, pred_labels)
    except:
        mcc = 0.0
    
    print(f"\n📈 VALIDATION METRICS:")
    print(f"   ✓ Accuracy:  {accuracy:.4f} ({accuracy*100:.2f}%)")
    print(f"   ✓ Precision: {precision:.4f}")
    print(f"   ✓ Recall:    {recall:.4f}")
    print(f"   ✓ F1-Score:  {f1:.4f}")
    print(f"   ✓ MCC:       {mcc:.4f}")
    
    # Per-residue accuracy
    correct_predictions = np.sum(true_labels == pred_labels)
    per_residue_acc = correct_predictions / len(true_labels)
    print(f"   ✓ Per-residue Accuracy: {per_residue_acc:.4f} ({per_residue_acc*100:.2f}%)")
    
    # Top-k accuracy
    top3_acc = np.mean([true_labels[i] in np.argsort(pred_probs_valid[i])[-3:] for i in range(len(true_labels))])
    top5_acc = np.mean([true_labels[i] in np.argsort(pred_probs_valid[i])[-5:] for i in range(len(true_labels))])
    
    print(f"   ✓ Top-3 Accuracy: {top3_acc:.4f} ({top3_acc*100:.2f}%)")
    print(f"   ✓ Top-5 Accuracy: {top5_acc:.4f} ({top5_acc*100:.2f}%)")
    
    # Confidence analysis
    max_probs = np.max(pred_probs_valid, axis=1)
    avg_confidence = np.mean(max_probs)
    print(f"   ✓ Average Confidence: {avg_confidence:.4f}")
    
    return {
        'accuracy': accuracy,
        'precision': precision,
        'recall': recall,
        'f1': f1,
        'mcc': mcc,
        'top3_accuracy': top3_acc,
        'top5_accuracy': top5_acc,
        'avg_confidence': avg_confidence,
        'true_labels': true_labels,
        'pred_labels': pred_labels,
        'pred_probs': pred_probs_valid,
        'sequence': aa_sequence,
        'total_residues': len(true_labels)
    }

def create_comprehensive_plots(results, protein_name):
    """Create comprehensive validation plots"""
    plt.style.use('default')
    
    # Create figure with multiple subplots
    fig = plt.figure(figsize=(20, 15))
    gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)
    
    fig.suptitle(f'GCNdesign Validation Analysis - {protein_name}', fontsize=16, fontweight='bold')
    
    # 1. Performance Metrics Bar Plot
    ax1 = fig.add_subplot(gs[0, 0])
    metrics = ['Accuracy', 'Precision', 'Recall', 'F1-Score', 'Top-3 Acc', 'Top-5 Acc']
    values = [results['accuracy'], results['precision'], results['recall'], 
              results['f1'], results['top3_accuracy'], results['top5_accuracy']]
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']
    
    bars = ax1.bar(metrics, values, color=colors, alpha=0.8)
    ax1.set_ylim(0, 1)
    ax1.set_ylabel('Score')
    ax1.set_title('Performance Metrics')
    ax1.grid(axis='y', alpha=0.3)
    plt.setp(ax1.get_xticklabels(), rotation=45, ha='right')
    
    for bar, value in zip(bars, values):
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01, 
                f'{value:.3f}', ha='center', va='bottom', fontsize=9)
    
    # 2. Confusion Matrix
    ax2 = fig.add_subplot(gs[0, 1])
    cm = confusion_matrix(results['true_labels'], results['pred_labels'])
    unique_labels = np.unique(np.concatenate([results['true_labels'], results['pred_labels']]))
    
    if len(unique_labels) <= 20:  # Only show if manageable size
        cm_labels = [amino_acids[i] for i in unique_labels]
        sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', ax=ax2,
                    xticklabels=cm_labels, yticklabels=cm_labels, cbar_kws={'shrink': 0.8})
        ax2.set_xlabel('Predicted')
        ax2.set_ylabel('True')
        ax2.set_title('Confusion Matrix')
    else:
        # For large matrices, show simplified version
        ax2.imshow(cm, cmap='Blues')
        ax2.set_title('Confusion Matrix (Simplified)')
        ax2.set_xlabel('Predicted Amino Acid Index')
        ax2.set_ylabel('True Amino Acid Index')
    
    # 3. Confidence Distribution
    ax3 = fig.add_subplot(gs[0, 2])
    max_probs = np.max(results['pred_probs'], axis=1)
    ax3.hist(max_probs, bins=30, alpha=0.7, color='skyblue', edgecolor='black')
    ax3.axvline(results['avg_confidence'], color='red', linestyle='--', 
                label=f'Mean: {results["avg_confidence"]:.3f}')
    ax3.set_xlabel('Prediction Confidence')
    ax3.set_ylabel('Frequency')
    ax3.set_title('Prediction Confidence Distribution')
    ax3.legend()
    ax3.grid(alpha=0.3)
    
    # 4. Per-Position Accuracy
    ax4 = fig.add_subplot(gs[1, :])
    correct_mask = (results['true_labels'] == results['pred_labels'])
    positions = range(1, len(correct_mask) + 1)
    colors = ['green' if correct else 'red' for correct in correct_mask]
    
    ax4.bar(positions, correct_mask.astype(int), color=colors, alpha=0.7, width=0.8)
    ax4.set_xlabel('Residue Position')
    ax4.set_ylabel('Correct (1) / Incorrect (0)')
    ax4.set_title('Per-Position Prediction Accuracy')
    ax4.set_ylim(-0.1, 1.1)
    ax4.grid(axis='y', alpha=0.3)
    
    # Add sequence on top
    sequence_str = ''.join([amino_acids[label] for label in results['true_labels']])
    if len(sequence_str) <= 100:  # Only show for reasonable lengths
        for i, aa in enumerate(sequence_str):
            ax4.text(i+1, 1.05, aa, ha='center', va='bottom', fontsize=8, 
                    color='green' if correct_mask[i] else 'red')
    
    # 5. Amino Acid Frequency Analysis
    ax5 = fig.add_subplot(gs[2, 0])
    true_counts = np.bincount(results['true_labels'], minlength=20)
    pred_counts = np.bincount(results['pred_labels'], minlength=20)
    
    x = np.arange(20)
    width = 0.35
    ax5.bar(x - width/2, true_counts, width, label='True', alpha=0.8, color='blue')
    ax5.bar(x + width/2, pred_counts, width, label='Predicted', alpha=0.8, color='orange')
    
    ax5.set_xlabel('Amino Acid')
    ax5.set_ylabel('Frequency')
    ax5.set_title('Amino Acid Distribution')
    ax5.set_xticks(x)
    ax5.set_xticklabels(amino_acids)
    ax5.legend()
    ax5.grid(axis='y', alpha=0.3)
    
    # 6. Model Performance Summary
    ax6 = fig.add_subplot(gs[2, 1])
    ax6.axis('off')
    
    summary_text = f"""
    MODEL PERFORMANCE SUMMARY
    
    Overall Accuracy: {results['accuracy']:.3f} ({results['accuracy']*100:.1f}%)
    
    Detailed Metrics:
    • Precision: {results['precision']:.3f}
    • Recall: {results['recall']:.3f}
    • F1-Score: {results['f1']:.3f}
    • MCC: {results['mcc']:.3f}
    
    Top-K Accuracy:
    • Top-3: {results['top3_accuracy']:.3f} ({results['top3_accuracy']*100:.1f}%)
    • Top-5: {results['top5_accuracy']:.3f} ({results['top5_accuracy']*100:.1f}%)
    
    Confidence:
    • Average: {results['avg_confidence']:.3f}
    
    Dataset:
    • Total Residues: {results['total_residues']}
    • Sequence Length: {len(results['sequence'])}
    """
    
    ax6.text(0.1, 0.5, summary_text, transform=ax6.transAxes, fontsize=10,
             verticalalignment='center', bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgray"))
    
    # 7. Prediction Quality by Amino Acid
    ax7 = fig.add_subplot(gs[2, 2])
    
    # Calculate per-amino-acid accuracy
    aa_accuracies = []
    aa_counts = []
    present_aas = []
    
    for aa_idx in range(20):
        mask = results['true_labels'] == aa_idx
        if np.sum(mask) > 0:  # Only include amino acids that appear
            acc = np.mean(results['pred_labels'][mask] == aa_idx)
            aa_accuracies.append(acc)
            aa_counts.append(np.sum(mask))
            present_aas.append(amino_acids[aa_idx])
    
    if aa_accuracies:
        bars = ax7.bar(range(len(present_aas)), aa_accuracies, 
                      color=plt.cm.viridis(np.array(aa_accuracies)))
        ax7.set_xlabel('Amino Acid')
        ax7.set_ylabel('Accuracy')
        ax7.set_title('Per-Amino Acid Prediction Accuracy')
        ax7.set_xticks(range(len(present_aas)))
        ax7.set_xticklabels(present_aas)
        ax7.set_ylim(0, 1)
        ax7.grid(axis='y', alpha=0.3)
        
        # Add count labels
        for i, (bar, count) in enumerate(zip(bars, aa_counts)):
            ax7.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02, 
                    f'n={count}', ha='center', va='bottom', fontsize=8)
    
    plt.tight_layout()
    plt.savefig(f'gcndesign_validation_{protein_name}.png', dpi=300, bbox_inches='tight')
    print(f"📊 Comprehensive plots saved as 'gcndesign_validation_{protein_name}.png'")
    
    return fig

def main():
    """Main analysis function"""
    print("🧬 GCNdesign Performance Analysis")
    print("=" * 50)
    
    # Use downloaded protein
    pdb_file = '1ubq.pdb'
    
    if not os.path.exists(pdb_file):
        print(f"❌ File {pdb_file} not found. Please download a protein structure.")
        return
    
    try:
        # Run analysis
        results = analyze_protein_predictions(pdb_file)
        
        # Create plots
        protein_name = os.path.splitext(pdb_file)[0]
        create_comprehensive_plots(results, protein_name)
        
        # Save results
        with open(f'validation_results_{protein_name}.pkl', 'wb') as f:
            pickle.dump(results, f)
        
        # Create summary report
        report = f"""
GCNdesign Validation Report - {protein_name.upper()}
=====================================================

PERFORMANCE METRICS:
• Overall Accuracy: {results['accuracy']:.4f} ({results['accuracy']*100:.2f}%)
• Precision: {results['precision']:.4f}
• Recall: {results['recall']:.4f}
• F1-Score: {results['f1']:.4f}
• Matthews Correlation Coefficient: {results['mcc']:.4f}

TOP-K ACCURACY:
• Top-1 (Standard): {results['accuracy']:.4f} ({results['accuracy']*100:.2f}%)
• Top-3: {results['top3_accuracy']:.4f} ({results['top3_accuracy']*100:.2f}%)
• Top-5: {results['top5_accuracy']:.4f} ({results['top5_accuracy']*100:.2f}%)

MODEL CONFIDENCE:
• Average Prediction Confidence: {results['avg_confidence']:.4f}

DATASET STATISTICS:
• Total Residues: {results['total_residues']}
• Protein Sequence: {' '.join(results['sequence'])}

FILES GENERATED:
• gcndesign_validation_{protein_name}.png - Comprehensive visualization
• validation_results_{protein_name}.pkl - Raw results data
• validation_report_{protein_name}.txt - This report

INTERPRETATION:
The model shows {'excellent' if results['accuracy'] > 0.8 else 'good' if results['accuracy'] > 0.6 else 'moderate'} performance with {results['accuracy']*100:.1f}% accuracy.
Top-3 accuracy of {results['top3_accuracy']*100:.1f}% indicates the model often includes the correct 
amino acid among its top 3 predictions, which is valuable for protein design applications.
"""
        
        with open(f'validation_report_{protein_name}.txt', 'w') as f:
            f.write(report)
        
        print(f"\n✅ Analysis complete! Generated files:")
        print(f"   📊 gcndesign_validation_{protein_name}.png")
        print(f"   📁 validation_results_{protein_name}.pkl") 
        print(f"   📄 validation_report_{protein_name}.txt")
        
        print(f"\n🎯 KEY RESULTS:")
        print(f"   Accuracy: {results['accuracy']*100:.1f}%")
        print(f"   Top-3 Accuracy: {results['top3_accuracy']*100:.1f}%")
        print(f"   Average Confidence: {results['avg_confidence']:.3f}")
        
    except Exception as e:
        print(f"❌ Error during analysis: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
