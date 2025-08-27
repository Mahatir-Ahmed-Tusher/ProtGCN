#!/usr/bin/env python3
"""
ProtGCN Overall Validation Metrics Analysis
Comprehensive evaluation of the trained model's performance across multiple proteins
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
    matthews_corrcoef, confusion_matrix, classification_report,
    top_k_accuracy_score
)
from tqdm import tqdm
import pickle
from datetime import datetime
import glob

# Add gcndesign to path
sys.path.append('.')
from gcndesign.predictor import Predictor
from gcndesign.hypara import HyperParam
from gcndesign.dataset import pdb2input

# Amino acid mapping
amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
               'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

class ProtGCNValidator:
    def __init__(self, device='cpu'):
        """Initialize the validator with the trained model"""
        self.device = device
        self.predictor = Predictor(device=device)
        self.hypara = HyperParam()
        
        # Storage for results
        self.results = []
        self.overall_metrics = {}
        
    def evaluate_single_protein(self, pdb_file):
        """Evaluate model performance on a single protein"""
        try:
            # Get predictions
            predictions = self.predictor.predict_logit_tensor(pdb_file)
            pred_probs = torch.softmax(torch.tensor(predictions), dim=1).numpy()
            
            # Get ground truth
            node, edgemat, adjmat, labels, mask, aa_sequence = pdb2input(pdb_file, self.hypara)
            
            # Filter valid residues
            valid_mask = mask.numpy()
            true_labels = labels.numpy()[valid_mask]
            pred_labels = np.argmax(predictions, axis=1)[valid_mask]
            pred_probs_valid = pred_probs[valid_mask]
            
            # Calculate metrics
            accuracy = accuracy_score(true_labels, pred_labels)
            precision = precision_score(true_labels, pred_labels, average='macro', zero_division=0)
            recall = recall_score(true_labels, pred_labels, average='macro', zero_division=0)
            f1 = f1_score(true_labels, pred_labels, average='macro', zero_division=0)
            
            try:
                mcc = matthews_corrcoef(true_labels, pred_labels)
            except:
                mcc = 0.0
            
            # Top-k accuracy (handle class mismatch)
            try:
                top3_acc = top_k_accuracy_score(true_labels, pred_probs_valid, k=3, labels=range(20))
                top5_acc = top_k_accuracy_score(true_labels, pred_probs_valid, k=5, labels=range(20))
            except:
                # Fallback: calculate manually
                top3_correct = 0
                top5_correct = 0
                for i, true_label in enumerate(true_labels):
                    top3_preds = np.argsort(pred_probs_valid[i])[-3:]
                    top5_preds = np.argsort(pred_probs_valid[i])[-5:]
                    if true_label in top3_preds:
                        top3_correct += 1
                    if true_label in top5_preds:
                        top5_correct += 1
                top3_acc = top3_correct / len(true_labels)
                top5_acc = top5_correct / len(true_labels)
            
            # Confidence analysis
            max_probs = np.max(pred_probs_valid, axis=1)
            avg_confidence = np.mean(max_probs)
            
            # Per-residue analysis
            correct_predictions = np.sum(true_labels == pred_labels)
            per_residue_acc = correct_predictions / len(true_labels)
            
            # Amino acid composition analysis
            true_seq = ''.join([amino_acids[i] for i in true_labels])
            pred_seq = ''.join([amino_acids[i] for i in pred_labels])
            
            result = {
                'protein_name': os.path.splitext(os.path.basename(pdb_file))[0],
                'pdb_file': pdb_file,
                'total_residues': len(true_labels),
                'correct_predictions': correct_predictions,
                'accuracy': accuracy,
                'precision': precision,
                'recall': recall,
                'f1_score': f1,
                'mcc': mcc,
                'top3_accuracy': top3_acc,
                'top5_accuracy': top5_acc,
                'avg_confidence': avg_confidence,
                'per_residue_accuracy': per_residue_acc,
                'true_sequence': true_seq,
                'predicted_sequence': pred_seq,
                'true_labels': true_labels,
                'pred_labels': pred_labels,
                'pred_probs': pred_probs_valid
            }
            
            return result
            
        except Exception as e:
            print(f"Error evaluating {pdb_file}: {e}")
            return None
    
    def evaluate_multiple_proteins(self, pdb_files, show_progress=True):
        """Evaluate model on multiple proteins"""
        print(f"🧬 Evaluating GCNdesign model on {len(pdb_files)} proteins...")
        
        if show_progress:
            pdb_files = tqdm(pdb_files, desc="Processing proteins")
        
        for pdb_file in pdb_files:
            result = self.evaluate_single_protein(pdb_file)
            if result:
                self.results.append(result)
        
        print(f"✅ Completed evaluation of {len(self.results)} proteins")
        return self.results
    
    def calculate_overall_metrics(self):
        """Calculate overall metrics across all evaluated proteins"""
        if not self.results:
            print("No results available. Run evaluation first.")
            return None
        
        # Aggregate metrics
        total_residues = sum(r['total_residues'] for r in self.results)
        total_correct = sum(r['correct_predictions'] for r in self.results)
        
        # Weighted averages (by number of residues)
        weighted_accuracy = sum(r['accuracy'] * r['total_residues'] for r in self.results) / total_residues
        weighted_precision = sum(r['precision'] * r['total_residues'] for r in self.results) / total_residues
        weighted_recall = sum(r['recall'] * r['total_residues'] for r in self.results) / total_residues
        weighted_f1 = sum(r['f1_score'] * r['total_residues'] for r in self.results) / total_residues
        weighted_mcc = sum(r['mcc'] * r['total_residues'] for r in self.results) / total_residues
        weighted_top3 = sum(r['top3_accuracy'] * r['total_residues'] for r in self.results) / total_residues
        weighted_top5 = sum(r['top5_accuracy'] * r['total_residues'] for r in self.results) / total_residues
        weighted_confidence = sum(r['avg_confidence'] * r['total_residues'] for r in self.results) / total_residues
        
        # Simple averages
        simple_accuracy = np.mean([r['accuracy'] for r in self.results])
        simple_precision = np.mean([r['precision'] for r in self.results])
        simple_recall = np.mean([r['recall'] for r in self.results])
        simple_f1 = np.mean([r['f1_score'] for r in self.results])
        
        # Overall per-residue accuracy
        overall_per_residue_acc = total_correct / total_residues
        
        # Collect all predictions for confusion matrix
        all_true_labels = np.concatenate([r['true_labels'] for r in self.results])
        all_pred_labels = np.concatenate([r['pred_labels'] for r in self.results])
        
        # Calculate overall confusion matrix
        overall_confusion = confusion_matrix(all_true_labels, all_pred_labels)
        
        self.overall_metrics = {
            'total_proteins': len(self.results),
            'total_residues': total_residues,
            'total_correct_predictions': total_correct,
            'overall_per_residue_accuracy': overall_per_residue_acc,
            
            # Weighted averages (more accurate for overall performance)
            'weighted_accuracy': weighted_accuracy,
            'weighted_precision': weighted_precision,
            'weighted_recall': weighted_recall,
            'weighted_f1_score': weighted_f1,
            'weighted_mcc': weighted_mcc,
            'weighted_top3_accuracy': weighted_top3,
            'weighted_top5_accuracy': weighted_top5,
            'weighted_avg_confidence': weighted_confidence,
            
            # Simple averages
            'simple_accuracy': simple_accuracy,
            'simple_precision': simple_precision,
            'simple_recall': simple_recall,
            'simple_f1_score': simple_f1,
            
            # Per-protein statistics
            'min_accuracy': min(r['accuracy'] for r in self.results),
            'max_accuracy': max(r['accuracy'] for r in self.results),
            'std_accuracy': np.std([r['accuracy'] for r in self.results]),
            
            # Confusion matrix
            'confusion_matrix': overall_confusion,
            'all_true_labels': all_true_labels,
            'all_pred_labels': all_pred_labels
        }
        
        return self.overall_metrics
    
    def print_overall_results(self):
        """Print comprehensive overall results"""
        if not self.overall_metrics:
            print("No overall metrics available. Run calculate_overall_metrics() first.")
            return
        
        print("\n" + "="*80)
        print("🧬 GCNdesign OVERALL VALIDATION METRICS")
        print("="*80)
        
        print(f"\n📊 DATASET SUMMARY:")
        print(f"   Total proteins evaluated: {self.overall_metrics['total_proteins']}")
        print(f"   Total residues analyzed: {self.overall_metrics['total_residues']:,}")
        print(f"   Total correct predictions: {self.overall_metrics['total_correct_predictions']:,}")
        
        print(f"\n🎯 OVERALL PERFORMANCE (Weighted by residue count):")
        print(f"   ✓ Per-residue Accuracy: {self.overall_metrics['overall_per_residue_accuracy']:.4f} ({self.overall_metrics['overall_per_residue_accuracy']*100:.2f}%)")
        print(f"   ✓ Weighted Accuracy:    {self.overall_metrics['weighted_accuracy']:.4f} ({self.overall_metrics['weighted_accuracy']*100:.2f}%)")
        print(f"   ✓ Weighted Precision:   {self.overall_metrics['weighted_precision']:.4f}")
        print(f"   ✓ Weighted Recall:      {self.overall_metrics['weighted_recall']:.4f}")
        print(f"   ✓ Weighted F1-Score:    {self.overall_metrics['weighted_f1_score']:.4f}")
        print(f"   ✓ Weighted MCC:         {self.overall_metrics['weighted_mcc']:.4f}")
        
        print(f"\n🏆 TOP-K ACCURACY:")
        print(f"   ✓ Top-3 Accuracy: {self.overall_metrics['weighted_top3_accuracy']:.4f} ({self.overall_metrics['weighted_top3_accuracy']*100:.2f}%)")
        print(f"   ✓ Top-5 Accuracy: {self.overall_metrics['weighted_top5_accuracy']:.4f} ({self.overall_metrics['weighted_top5_accuracy']*100:.2f}%)")
        
        print(f"\n🎲 CONFIDENCE ANALYSIS:")
        print(f"   ✓ Average Confidence: {self.overall_metrics['weighted_avg_confidence']:.4f}")
        
        print(f"\n📈 PER-PROTEIN STATISTICS:")
        print(f"   ✓ Min Accuracy: {self.overall_metrics['min_accuracy']:.4f} ({self.overall_metrics['min_accuracy']*100:.2f}%)")
        print(f"   ✓ Max Accuracy: {self.overall_metrics['max_accuracy']:.4f} ({self.overall_metrics['max_accuracy']*100:.2f}%)")
        print(f"   ✓ Std Accuracy: {self.overall_metrics['std_accuracy']:.4f} ({self.overall_metrics['std_accuracy']*100:.2f}%)")
        
        # Performance interpretation
        accuracy = self.overall_metrics['weighted_accuracy']
        if accuracy > 0.8:
            performance_level = "EXCELLENT"
        elif accuracy > 0.6:
            performance_level = "GOOD"
        elif accuracy > 0.4:
            performance_level = "MODERATE"
        else:
            performance_level = "NEEDS IMPROVEMENT"
        
        print(f"\n🏅 PERFORMANCE ASSESSMENT:")
        print(f"   Model Performance Level: {performance_level}")
        print(f"   Overall Accuracy: {accuracy*100:.1f}%")
        
        if self.overall_metrics['weighted_top3_accuracy'] > 0.7:
            print(f"   ✓ Top-3 accuracy > 70% - Good for protein design applications")
        
        if self.overall_metrics['weighted_avg_confidence'] > 0.5:
            print(f"   ✓ Average confidence > 50% - Model shows reasonable certainty")
    
    def create_comprehensive_plots(self, save_dir="validation_results"):
        """Create comprehensive validation plots"""
        os.makedirs(save_dir, exist_ok=True)
        
        # Create figure with multiple subplots
        fig = plt.figure(figsize=(20, 16))
        
        # 1. Overall accuracy distribution
        plt.subplot(3, 3, 1)
        accuracies = [r['accuracy'] for r in self.results]
        plt.hist(accuracies, bins=20, alpha=0.7, color='skyblue', edgecolor='black')
        plt.axvline(np.mean(accuracies), color='red', linestyle='--', label=f'Mean: {np.mean(accuracies):.3f}')
        plt.xlabel('Accuracy')
        plt.ylabel('Number of Proteins')
        plt.title('Accuracy Distribution Across Proteins')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # 2. Per-residue accuracy vs protein size
        plt.subplot(3, 3, 2)
        sizes = [r['total_residues'] for r in self.results]
        plt.scatter(sizes, accuracies, alpha=0.6, color='green')
        plt.xlabel('Protein Size (residues)')
        plt.ylabel('Accuracy')
        plt.title('Accuracy vs Protein Size')
        plt.grid(True, alpha=0.3)
        
        # 3. Top-k accuracy comparison
        plt.subplot(3, 3, 3)
        top1_acc = [r['accuracy'] for r in self.results]
        top3_acc = [r['top3_accuracy'] for r in self.results]
        top5_acc = [r['top5_accuracy'] for r in self.results]
        
        x = np.arange(len(self.results))
        width = 0.25
        plt.bar(x - width, top1_acc, width, label='Top-1', alpha=0.8)
        plt.bar(x, top3_acc, width, label='Top-3', alpha=0.8)
        plt.bar(x + width, top5_acc, width, label='Top-5', alpha=0.8)
        plt.xlabel('Proteins')
        plt.ylabel('Accuracy')
        plt.title('Top-K Accuracy Comparison')
        plt.legend()
        plt.xticks(x[::5], [r['protein_name'] for r in self.results[::5]], rotation=45)
        
        # 4. Confusion matrix
        plt.subplot(3, 3, 4)
        cm = self.overall_metrics['confusion_matrix']
        sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', 
                   xticklabels=amino_acids, yticklabels=amino_acids)
        plt.title('Overall Confusion Matrix')
        plt.xlabel('Predicted')
        plt.ylabel('True')
        
        # 5. Confidence distribution
        plt.subplot(3, 3, 5)
        confidences = [r['avg_confidence'] for r in self.results]
        plt.hist(confidences, bins=20, alpha=0.7, color='orange', edgecolor='black')
        plt.axvline(np.mean(confidences), color='red', linestyle='--', label=f'Mean: {np.mean(confidences):.3f}')
        plt.xlabel('Average Confidence')
        plt.ylabel('Number of Proteins')
        plt.title('Confidence Distribution')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # 6. Metrics comparison
        plt.subplot(3, 3, 6)
        metrics = ['Accuracy', 'Precision', 'Recall', 'F1-Score']
        weighted_values = [
            self.overall_metrics['weighted_accuracy'],
            self.overall_metrics['weighted_precision'],
            self.overall_metrics['weighted_recall'],
            self.overall_metrics['weighted_f1_score']
        ]
        simple_values = [
            self.overall_metrics['simple_accuracy'],
            self.overall_metrics['simple_precision'],
            self.overall_metrics['simple_recall'],
            self.overall_metrics['simple_f1_score']
        ]
        
        x = np.arange(len(metrics))
        width = 0.35
        plt.bar(x - width/2, weighted_values, width, label='Weighted', alpha=0.8)
        plt.bar(x + width/2, simple_values, width, label='Simple Average', alpha=0.8)
        plt.xlabel('Metrics')
        plt.ylabel('Score')
        plt.title('Overall Metrics Comparison')
        plt.xticks(x, metrics)
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # 7. Protein-wise performance ranking
        plt.subplot(3, 3, 7)
        sorted_results = sorted(self.results, key=lambda x: x['accuracy'], reverse=True)
        top_proteins = sorted_results[:10]
        protein_names = [r['protein_name'] for r in top_proteins]
        protein_accuracies = [r['accuracy'] for r in top_proteins]
        
        plt.barh(range(len(protein_names)), protein_accuracies, color='lightcoral')
        plt.yticks(range(len(protein_names)), protein_names)
        plt.xlabel('Accuracy')
        plt.title('Top 10 Proteins by Accuracy')
        plt.grid(True, alpha=0.3)
        
        # 8. Sequence length vs performance
        plt.subplot(3, 3, 8)
        plt.scatter(sizes, confidences, alpha=0.6, color='purple')
        plt.xlabel('Protein Size (residues)')
        plt.ylabel('Average Confidence')
        plt.title('Confidence vs Protein Size')
        plt.grid(True, alpha=0.3)
        
        # 9. Summary statistics
        plt.subplot(3, 3, 9)
        summary_stats = [
            f"Total Proteins: {self.overall_metrics['total_proteins']}",
            f"Total Residues: {self.overall_metrics['total_residues']:,}",
            f"Overall Accuracy: {self.overall_metrics['weighted_accuracy']*100:.1f}%",
            f"Top-3 Accuracy: {self.overall_metrics['weighted_top3_accuracy']*100:.1f}%",
            f"Avg Confidence: {self.overall_metrics['weighted_avg_confidence']:.3f}",
            f"Min Accuracy: {self.overall_metrics['min_accuracy']*100:.1f}%",
            f"Max Accuracy: {self.overall_metrics['max_accuracy']*100:.1f}%"
        ]
        
        plt.text(0.1, 0.9, '\n'.join(summary_stats), transform=plt.gca().transAxes,
                fontsize=10, verticalalignment='top', fontfamily='monospace',
                bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgray", alpha=0.8))
        plt.axis('off')
        plt.title('Summary Statistics')
        
        plt.tight_layout()
        plt.savefig(os.path.join(save_dir, 'overall_validation_metrics.png'), 
                   dpi=300, bbox_inches='tight')
        plt.show()
        
        print(f"📊 Comprehensive plots saved to: {save_dir}/overall_validation_metrics.png")
    
    def save_results(self, save_dir="validation_results"):
        """Save all results to files"""
        os.makedirs(save_dir, exist_ok=True)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # Save detailed results
        results_df = pd.DataFrame(self.results)
        results_df.to_csv(os.path.join(save_dir, f'detailed_results_{timestamp}.csv'), index=False)
        
        # Save overall metrics
        with open(os.path.join(save_dir, f'overall_metrics_{timestamp}.pkl'), 'wb') as f:
            pickle.dump(self.overall_metrics, f)
        
        # Create comprehensive report
        report = f"""
GCNdesign Overall Validation Report
==================================
Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}

DATASET SUMMARY:
- Total proteins evaluated: {self.overall_metrics['total_proteins']}
- Total residues analyzed: {self.overall_metrics['total_residues']:,}
- Total correct predictions: {self.overall_metrics['total_correct_predictions']:,}

OVERALL PERFORMANCE (Weighted by residue count):
- Per-residue Accuracy: {self.overall_metrics['overall_per_residue_accuracy']:.4f} ({self.overall_metrics['overall_per_residue_accuracy']*100:.2f}%)
- Weighted Accuracy: {self.overall_metrics['weighted_accuracy']:.4f} ({self.overall_metrics['weighted_accuracy']*100:.2f}%)
- Weighted Precision: {self.overall_metrics['weighted_precision']:.4f}
- Weighted Recall: {self.overall_metrics['weighted_recall']:.4f}
- Weighted F1-Score: {self.overall_metrics['weighted_f1_score']:.4f}
- Weighted MCC: {self.overall_metrics['weighted_mcc']:.4f}

TOP-K ACCURACY:
- Top-3 Accuracy: {self.overall_metrics['weighted_top3_accuracy']:.4f} ({self.overall_metrics['weighted_top3_accuracy']*100:.2f}%)
- Top-5 Accuracy: {self.overall_metrics['weighted_top5_accuracy']:.4f} ({self.overall_metrics['weighted_top5_accuracy']*100:.2f}%)

CONFIDENCE ANALYSIS:
- Average Confidence: {self.overall_metrics['weighted_avg_confidence']:.4f}

PER-PROTEIN STATISTICS:
- Min Accuracy: {self.overall_metrics['min_accuracy']:.4f} ({self.overall_metrics['min_accuracy']*100:.2f}%)
- Max Accuracy: {self.overall_metrics['max_accuracy']:.4f} ({self.overall_metrics['max_accuracy']*100:.2f}%)
- Std Accuracy: {self.overall_metrics['std_accuracy']:.4f} ({self.overall_metrics['std_accuracy']*100:.2f}%)

PERFORMANCE ASSESSMENT:
- Model Performance Level: {'EXCELLENT' if self.overall_metrics['weighted_accuracy'] > 0.8 else 'GOOD' if self.overall_metrics['weighted_accuracy'] > 0.6 else 'MODERATE' if self.overall_metrics['weighted_accuracy'] > 0.4 else 'NEEDS IMPROVEMENT'}
- Overall Accuracy: {self.overall_metrics['weighted_accuracy']*100:.1f}%

GENERATED FILES:
- detailed_results_{timestamp}.csv - Detailed per-protein results
- overall_metrics_{timestamp}.pkl - Overall metrics data
- overall_validation_metrics.png - Comprehensive visualization plots
- validation_report_{timestamp}.txt - This report
"""
        
        with open(os.path.join(save_dir, f'validation_report_{timestamp}.txt'), 'w') as f:
            f.write(report)
        
        print(f"💾 Results saved to: {save_dir}/")
        print(f"   📊 detailed_results_{timestamp}.csv")
        print(f"   📁 overall_metrics_{timestamp}.pkl")
        print(f"   📄 validation_report_{timestamp}.txt")

def main():
    """Main function to run comprehensive validation"""
    print("🧬 ProtGCN Overall Validation Metrics Analysis")
    print("=" * 60)
    
    # Initialize validator
    validator = ProtGCNValidator(device='cpu')
    
    # Find PDB files to evaluate
    # You can modify this to use your own test dataset
    pdb_files = []
    
    # Option 1: Use built-in example proteins
    example_proteins = ['1ubq.pdb', '1ZNI.pdb', '1AKI.pdb']
    for protein in example_proteins:
        if os.path.exists(protein):
            pdb_files.append(protein)
    
    # Option 2: Use all PDB files in current directory
    if not pdb_files:
        pdb_files = glob.glob("*.pdb")
    
    # Option 3: Use specific test dataset
    if not pdb_files:
        print("No PDB files found. Please provide a list of PDB files to evaluate.")
        print("You can:")
        print("1. Place PDB files in the current directory")
        print("2. Modify the script to point to your test dataset")
        print("3. Use the built-in example proteins")
        return
    
    print(f"📁 Found {len(pdb_files)} PDB files to evaluate")
    
    # Evaluate proteins
    results = validator.evaluate_multiple_proteins(pdb_files)
    
    if not results:
        print("❌ No proteins were successfully evaluated")
        return
    
    # Calculate overall metrics
    overall_metrics = validator.calculate_overall_metrics()
    
    # Print results
    validator.print_overall_results()
    
    # Create plots
    print("\n📊 Creating comprehensive plots...")
    validator.create_comprehensive_plots()
    
    # Save results
    print("\n💾 Saving results...")
    validator.save_results()
    
    print("\n✅ Overall validation analysis complete!")

if __name__ == "__main__":
    main()
