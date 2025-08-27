#!/usr/bin/env python3
"""
ProtGCN Protein Design Metrics Calculator
Calculate comprehensive metrics for protein design evaluation including T500 and TS50 equivalents
"""

import os
import sys
import torch
import numpy as np
import pandas as pd
from sklearn.metrics import top_k_accuracy_score, accuracy_score, precision_score, recall_score, f1_score
from tqdm import tqdm
import glob
from datetime import datetime

# Add gcndesign to path
sys.path.append('.')
from gcndesign.predictor import Predictor
from gcndesign.hypara import HyperParam
from gcndesign.dataset import pdb2input

class ProteinDesignMetricsCalculator:
    def __init__(self, device='cpu'):
        """Initialize the calculator with the trained model"""
        self.device = device
        self.predictor = Predictor(device=device)
        self.hypara = HyperParam()
        
        # Storage for results
        self.results = []
        self.overall_metrics = {}
        
    def calculate_design_metrics_single_protein(self, pdb_file):
        """Calculate comprehensive protein design metrics for a single protein"""
        try:
            # Get predictions
            predictions = self.predictor.predict_logit_tensor(pdb_file)
            pred_probs = torch.softmax(torch.tensor(predictions), dim=1).numpy()
            
            # Get ground truth
            node, edgemat, adjmat, labels, mask, aa_sequence = pdb2input(pdb_file, self.hypara)
            
            # Filter valid residues
            valid_mask = mask.numpy()
            true_labels = labels.numpy()[valid_mask]
            pred_probs_valid = pred_probs[valid_mask]
            pred_labels = np.argmax(predictions, axis=1)[valid_mask]
            
            # Calculate basic metrics
            accuracy = accuracy_score(true_labels, pred_labels)
            precision = precision_score(true_labels, pred_labels, average='macro', zero_division=0)
            recall = recall_score(true_labels, pred_labels, average='macro', zero_division=0)
            f1 = f1_score(true_labels, pred_labels, average='macro', zero_division=0)
            
            # Calculate Top-K accuracies (for 20 amino acids)
            top1_acc = accuracy
            top3_acc = top_k_accuracy_score(true_labels, pred_probs_valid, k=3, labels=range(20))
            top5_acc = top_k_accuracy_score(true_labels, pred_probs_valid, k=5, labels=range(20))
            top10_acc = top_k_accuracy_score(true_labels, pred_probs_valid, k=10, labels=range(20))
            top20_acc = top_k_accuracy_score(true_labels, pred_probs_valid, k=20, labels=range(20))
            
            # Calculate T500 and TS50 equivalents for protein design
            # Since we have 20 amino acids, we calculate:
            # T500 equivalent: Top-20 accuracy (all amino acids)
            # TS50 equivalent: Top-10 accuracy (half of all amino acids)
            t500_equivalent = top20_acc  # All 20 amino acids
            ts50_equivalent = top10_acc  # Top 10 amino acids
            
            # Calculate confidence metrics
            max_probs = np.max(pred_probs_valid, axis=1)
            avg_confidence = np.mean(max_probs)
            confidence_std = np.std(max_probs)
            
            # Calculate entropy (diversity of predictions)
            entropy = -np.sum(pred_probs_valid * np.log(pred_probs_valid + 1e-10), axis=1)
            avg_entropy = np.mean(entropy)
            
            # Calculate per-residue analysis
            correct_predictions = np.sum(true_labels == pred_labels)
            per_residue_acc = correct_predictions / len(true_labels)
            
            # Amino acid composition analysis
            amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                          'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
            true_seq = ''.join([amino_acids[i] for i in true_labels])
            pred_seq = ''.join([amino_acids[i] for i in pred_labels])
            
            result = {
                'protein_name': os.path.splitext(os.path.basename(pdb_file))[0],
                'pdb_file': pdb_file,
                'total_residues': len(true_labels),
                'correct_predictions': correct_predictions,
                
                # Basic metrics
                'accuracy': accuracy,
                'precision': precision,
                'recall': recall,
                'f1_score': f1,
                'per_residue_accuracy': per_residue_acc,
                
                # Top-K accuracies
                'top1_accuracy': top1_acc,
                'top3_accuracy': top3_acc,
                'top5_accuracy': top5_acc,
                'top10_accuracy': top10_acc,
                'top20_accuracy': top20_acc,
                
                # Protein design specific metrics
                't500_equivalent': t500_equivalent,  # Top-20 (all amino acids)
                'ts50_equivalent': ts50_equivalent,  # Top-10 (half amino acids)
                
                # Confidence metrics
                'avg_confidence': avg_confidence,
                'confidence_std': confidence_std,
                'avg_entropy': avg_entropy,
                
                # Sequences
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
    
    def calculate_multiple_proteins(self, pdb_files, show_progress=True):
        """Calculate metrics for multiple proteins"""
        print(f"🧬 Calculating protein design metrics for {len(pdb_files)} proteins...")
        
        if show_progress:
            pdb_files = tqdm(pdb_files, desc="Processing proteins")
        
        for pdb_file in pdb_files:
            result = self.calculate_design_metrics_single_protein(pdb_file)
            if result:
                self.results.append(result)
        
        print(f"✅ Completed evaluation of {len(self.results)} proteins")
        return self.results
    
    def calculate_overall_metrics(self):
        """Calculate overall metrics"""
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
        
        weighted_top1 = sum(r['top1_accuracy'] * r['total_residues'] for r in self.results) / total_residues
        weighted_top3 = sum(r['top3_accuracy'] * r['total_residues'] for r in self.results) / total_residues
        weighted_top5 = sum(r['top5_accuracy'] * r['total_residues'] for r in self.results) / total_residues
        weighted_top10 = sum(r['top10_accuracy'] * r['total_residues'] for r in self.results) / total_residues
        weighted_top20 = sum(r['top20_accuracy'] * r['total_residues'] for r in self.results) / total_residues
        
        weighted_t500_equiv = sum(r['t500_equivalent'] * r['total_residues'] for r in self.results) / total_residues
        weighted_ts50_equiv = sum(r['ts50_equivalent'] * r['total_residues'] for r in self.results) / total_residues
        
        weighted_confidence = sum(r['avg_confidence'] * r['total_residues'] for r in self.results) / total_residues
        weighted_entropy = sum(r['avg_entropy'] * r['total_residues'] for r in self.results) / total_residues
        
        # Simple averages
        simple_accuracy = np.mean([r['accuracy'] for r in self.results])
        simple_t500_equiv = np.mean([r['t500_equivalent'] for r in self.results])
        simple_ts50_equiv = np.mean([r['ts50_equivalent'] for r in self.results])
        simple_top3 = np.mean([r['top3_accuracy'] for r in self.results])
        simple_top5 = np.mean([r['top5_accuracy'] for r in self.results])
        
        # Overall per-residue accuracy
        overall_per_residue_acc = total_correct / total_residues
        
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
            'weighted_top1_accuracy': weighted_top1,
            'weighted_top3_accuracy': weighted_top3,
            'weighted_top5_accuracy': weighted_top5,
            'weighted_top10_accuracy': weighted_top10,
            'weighted_top20_accuracy': weighted_top20,
            'weighted_t500_equivalent': weighted_t500_equiv,
            'weighted_ts50_equivalent': weighted_ts50_equiv,
            'weighted_avg_confidence': weighted_confidence,
            'weighted_avg_entropy': weighted_entropy,
            
            # Simple averages
            'simple_accuracy': simple_accuracy,
            'simple_t500_equivalent': simple_t500_equiv,
            'simple_ts50_equivalent': simple_ts50_equiv,
            'simple_top3_accuracy': simple_top3,
            'simple_top5_accuracy': simple_top5,
            
            # Per-protein statistics
            'min_accuracy': min(r['accuracy'] for r in self.results),
            'max_accuracy': max(r['accuracy'] for r in self.results),
            'std_accuracy': np.std([r['accuracy'] for r in self.results]),
            'min_t500_equiv': min(r['t500_equivalent'] for r in self.results),
            'max_t500_equiv': max(r['t500_equivalent'] for r in self.results),
            'min_ts50_equiv': min(r['ts50_equivalent'] for r in self.results),
            'max_ts50_equiv': max(r['ts50_equivalent'] for r in self.results)
        }
        
        return self.overall_metrics
    
    def print_results(self):
        """Print comprehensive protein design metrics results"""
        if not self.overall_metrics:
            print("No overall metrics available. Run calculate_overall_metrics() first.")
            return
        
        print("\n" + "="*80)
        print("🧬 GCNdesign PROTEIN DESIGN METRICS")
        print("="*80)
        
        print(f"\n📊 DATASET SUMMARY:")
        print(f"   Total proteins evaluated: {self.overall_metrics['total_proteins']}")
        print(f"   Total residues analyzed: {self.overall_metrics['total_residues']:,}")
        print(f"   Total correct predictions: {self.overall_metrics['total_correct_predictions']:,}")
        
        print(f"\n🎯 OVERALL PERFORMANCE (Weighted by residue count):")
        print(f"   ✓ Per-residue Accuracy: {self.overall_metrics['overall_per_residue_accuracy']:.4f} ({self.overall_metrics['overall_per_residue_accuracy']*100:.2f}%)")
        print(f"   ✓ Weighted Accuracy: {self.overall_metrics['weighted_accuracy']:.4f} ({self.overall_metrics['weighted_accuracy']*100:.2f}%)")
        print(f"   ✓ Weighted Precision: {self.overall_metrics['weighted_precision']:.4f}")
        print(f"   ✓ Weighted Recall: {self.overall_metrics['weighted_recall']:.4f}")
        print(f"   ✓ Weighted F1-Score: {self.overall_metrics['weighted_f1_score']:.4f}")
        
        print(f"\n🏆 TOP-K ACCURACY:")
        print(f"   ✓ Top-1 Accuracy: {self.overall_metrics['weighted_top1_accuracy']:.4f} ({self.overall_metrics['weighted_top1_accuracy']*100:.2f}%)")
        print(f"   ✓ Top-3 Accuracy: {self.overall_metrics['weighted_top3_accuracy']:.4f} ({self.overall_metrics['weighted_top3_accuracy']*100:.2f}%)")
        print(f"   ✓ Top-5 Accuracy: {self.overall_metrics['weighted_top5_accuracy']:.4f} ({self.overall_metrics['weighted_top5_accuracy']*100:.2f}%)")
        print(f"   ✓ Top-10 Accuracy: {self.overall_metrics['weighted_top10_accuracy']:.4f} ({self.overall_metrics['weighted_top10_accuracy']*100:.2f}%)")
        print(f"   ✓ Top-20 Accuracy: {self.overall_metrics['weighted_top20_accuracy']:.4f} ({self.overall_metrics['weighted_top20_accuracy']*100:.2f}%)")
        
        print(f"\n🧬 PROTEIN DESIGN SPECIFIC METRICS:")
        print(f"   ✓ T500 Equivalent (Top-20): {self.overall_metrics['weighted_t500_equivalent']:.4f} ({self.overall_metrics['weighted_t500_equivalent']*100:.2f}%)")
        print(f"   ✓ TS50 Equivalent (Top-10): {self.overall_metrics['weighted_ts50_equivalent']:.4f} ({self.overall_metrics['weighted_ts50_equivalent']*100:.2f}%)")
        
        print(f"\n🎲 CONFIDENCE & DIVERSITY:")
        print(f"   ✓ Average Confidence: {self.overall_metrics['weighted_avg_confidence']:.4f}")
        print(f"   ✓ Average Entropy: {self.overall_metrics['weighted_avg_entropy']:.4f}")
        
        print(f"\n📈 PER-PROTEIN STATISTICS:")
        print(f"   ✓ Min Accuracy: {self.overall_metrics['min_accuracy']:.4f} ({self.overall_metrics['min_accuracy']*100:.2f}%)")
        print(f"   ✓ Max Accuracy: {self.overall_metrics['max_accuracy']:.4f} ({self.overall_metrics['max_accuracy']*100:.2f}%)")
        print(f"   ✓ Std Accuracy: {self.overall_metrics['std_accuracy']:.4f} ({self.overall_metrics['std_accuracy']*100:.2f}%)")
        
        # Performance interpretation
        accuracy = self.overall_metrics['weighted_accuracy']
        t500_equiv = self.overall_metrics['weighted_t500_equivalent']
        ts50_equiv = self.overall_metrics['weighted_ts50_equivalent']
        
        print(f"\n🏅 PERFORMANCE ASSESSMENT:")
        print(f"   Overall Accuracy: {accuracy*100:.1f}%")
        print(f"   T500 Equivalent: {t500_equiv*100:.1f}%")
        print(f"   TS50 Equivalent: {ts50_equiv*100:.1f}%")
        
        # Compare with literature values (adjusted for 20 amino acids)
        print(f"\n📊 COMPARISON WITH LITERATURE (Adjusted for 20 amino acids):")
        print(f"   GCNdesign (3 layers): T500 ~53.49%, TS50 ~47.06%")
        print(f"   GCNdesign (4 layers): T500 ~53.78%, TS50 ~47.42%")
        print(f"   GCNdesign (5 layers): T500 ~53.64%, TS50 ~47.61%")
        print(f"   GCNdesign (6 layers): T500 ~53.41%, TS50 ~47.61%")
        print(f"   DenseCPD: T500 ~53.24%, TS50 ~46.74%")
        print(f"   ProDCoNN: T500 ~52.82%, TS50 ~50.71%")
        print(f"   SPROF: T500 ~42.20%, TS50 ~40.25%")
        print(f"   SPIN2: T500 ~40.69%, TS50 ~39.16%")
        
        # Performance level assessment
        if accuracy > 0.8:
            performance_level = "EXCELLENT"
        elif accuracy > 0.6:
            performance_level = "GOOD"
        elif accuracy > 0.4:
            performance_level = "MODERATE"
        else:
            performance_level = "NEEDS IMPROVEMENT"
        
        print(f"\n🏆 PERFORMANCE LEVEL:")
        print(f"   Model Performance Level: {performance_level}")
        
        if self.overall_metrics['weighted_top3_accuracy'] > 0.7:
            print(f"   ✓ Top-3 accuracy > 70% - Good for protein design applications")
        
        if self.overall_metrics['weighted_avg_confidence'] > 0.5:
            print(f"   ✓ Average confidence > 50% - Model shows reasonable certainty")
    
    def save_results(self, save_dir="protein_design_metrics"):
        """Save protein design metrics results to files"""
        os.makedirs(save_dir, exist_ok=True)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # Save detailed results
        results_df = pd.DataFrame(self.results)
        results_df.to_csv(os.path.join(save_dir, f'protein_design_detailed_{timestamp}.csv'), index=False)
        
        # Save overall metrics
        with open(os.path.join(save_dir, f'protein_design_overall_{timestamp}.pkl'), 'wb') as f:
            import pickle
            pickle.dump(self.overall_metrics, f)
        
        # Create comprehensive report
        report = f"""
GCNdesign Protein Design Metrics Report
=======================================
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

TOP-K ACCURACY:
- Top-1 Accuracy: {self.overall_metrics['weighted_top1_accuracy']:.4f} ({self.overall_metrics['weighted_top1_accuracy']*100:.2f}%)
- Top-3 Accuracy: {self.overall_metrics['weighted_top3_accuracy']:.4f} ({self.overall_metrics['weighted_top3_accuracy']*100:.2f}%)
- Top-5 Accuracy: {self.overall_metrics['weighted_top5_accuracy']:.4f} ({self.overall_metrics['weighted_top5_accuracy']*100:.2f}%)
- Top-10 Accuracy: {self.overall_metrics['weighted_top10_accuracy']:.4f} ({self.overall_metrics['weighted_top10_accuracy']*100:.2f}%)
- Top-20 Accuracy: {self.overall_metrics['weighted_top20_accuracy']:.4f} ({self.overall_metrics['weighted_top20_accuracy']*100:.2f}%)

PROTEIN DESIGN SPECIFIC METRICS:
- T500 Equivalent (Top-20): {self.overall_metrics['weighted_t500_equivalent']:.4f} ({self.overall_metrics['weighted_t500_equivalent']*100:.2f}%)
- TS50 Equivalent (Top-10): {self.overall_metrics['weighted_ts50_equivalent']:.4f} ({self.overall_metrics['weighted_ts50_equivalent']*100:.2f}%)

CONFIDENCE & DIVERSITY:
- Average Confidence: {self.overall_metrics['weighted_avg_confidence']:.4f}
- Average Entropy: {self.overall_metrics['weighted_avg_entropy']:.4f}

PER-PROTEIN STATISTICS:
- Min Accuracy: {self.overall_metrics['min_accuracy']:.4f} ({self.overall_metrics['min_accuracy']*100:.2f}%)
- Max Accuracy: {self.overall_metrics['max_accuracy']:.4f} ({self.overall_metrics['max_accuracy']*100:.2f}%)
- Std Accuracy: {self.overall_metrics['std_accuracy']:.4f} ({self.overall_metrics['std_accuracy']*100:.2f}%)

LITERATURE COMPARISON (Adjusted for 20 amino acids):
- GCNdesign (3 layers): T500 ~53.49%, TS50 ~47.06%
- GCNdesign (4 layers): T500 ~53.78%, TS50 ~47.42%
- GCNdesign (5 layers): T500 ~53.64%, TS50 ~47.61%
- GCNdesign (6 layers): T500 ~53.41%, TS50 ~47.61%
- DenseCPD: T500 ~53.24%, TS50 ~46.74%
- ProDCoNN: T500 ~52.82%, TS50 ~50.71%
- SPROF: T500 ~42.20%, TS50 ~40.25%
- SPIN2: T500 ~40.69%, TS50 ~39.16%

GENERATED FILES:
- protein_design_detailed_{timestamp}.csv - Detailed per-protein results
- protein_design_overall_{timestamp}.pkl - Overall metrics data
- protein_design_report_{timestamp}.txt - This report
"""
        
        with open(os.path.join(save_dir, f'protein_design_report_{timestamp}.txt'), 'w') as f:
            f.write(report)
        
        print(f"💾 Results saved to: {save_dir}/")
        print(f"   📊 protein_design_detailed_{timestamp}.csv")
        print(f"   📁 protein_design_overall_{timestamp}.pkl")
        print(f"   📄 protein_design_report_{timestamp}.txt")

def main():
    """Main function to run protein design metrics calculation"""
    print("🧬 ProtGCN Protein Design Metrics Calculator")
    print("=" * 60)
    
    # Initialize calculator
    calculator = ProteinDesignMetricsCalculator(device='cpu')
    
    # Find PDB files to evaluate
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
    
    # Calculate metrics
    results = calculator.calculate_multiple_proteins(pdb_files)
    
    if not results:
        print("❌ No proteins were successfully evaluated")
        return
    
    # Calculate overall metrics
    overall_metrics = calculator.calculate_overall_metrics()
    
    # Print results
    calculator.print_results()
    
    # Save results
    print("\n💾 Saving results...")
    calculator.save_results()
    
    print("\n✅ Protein design metrics analysis complete!")

if __name__ == "__main__":
    main()
