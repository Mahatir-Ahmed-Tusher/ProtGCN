#!/usr/bin/env python3
"""
ProtGCN T500 and TS50 Metrics Calculator
Calculate Top-500 and Top-50 accuracy metrics for protein design evaluation
"""

import os
import sys
import torch
import numpy as np
import pandas as pd
from sklearn.metrics import top_k_accuracy_score
from tqdm import tqdm
import glob
from datetime import datetime

# Add gcndesign to path
sys.path.append('.')
from gcndesign.predictor import Predictor
from gcndesign.hypara import HyperParam
from gcndesign.dataset import pdb2input

class T500TS50Calculator:
    def __init__(self, device='cpu'):
        """Initialize the calculator with the trained model"""
        self.device = device
        self.predictor = Predictor(device=device)
        self.hypara = HyperParam()
        
        # Storage for results
        self.results = []
        self.overall_metrics = {}
        
    def calculate_t500_ts50_single_protein(self, pdb_file):
        """Calculate T500 and TS50 metrics for a single protein"""
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
            
            # Calculate T500 and TS50
            try:
                t500_acc = top_k_accuracy_score(true_labels, pred_probs_valid, k=500, labels=range(20))
            except:
                # Fallback: calculate manually for T500
                t500_correct = 0
                for i, true_label in enumerate(true_labels):
                    top500_preds = np.argsort(pred_probs_valid[i])[-500:]
                    if true_label in top500_preds:
                        t500_correct += 1
                t500_acc = t500_correct / len(true_labels)
            
            try:
                ts50_acc = top_k_accuracy_score(true_labels, pred_probs_valid, k=50, labels=range(20))
            except:
                # Fallback: calculate manually for TS50
                ts50_correct = 0
                for i, true_label in enumerate(true_labels):
                    top50_preds = np.argsort(pred_probs_valid[i])[-50:]
                    if true_label in top50_preds:
                        ts50_correct += 1
                ts50_acc = ts50_correct / len(true_labels)
            
            # Also calculate other common metrics for comparison
            pred_labels = np.argmax(predictions, axis=1)[valid_mask]
            accuracy = np.mean(true_labels == pred_labels)
            
            # Top-1, Top-3, Top-5, Top-10, Top-20
            top1_acc = accuracy
            top3_acc = top_k_accuracy_score(true_labels, pred_probs_valid, k=3, labels=range(20))
            top5_acc = top_k_accuracy_score(true_labels, pred_probs_valid, k=5, labels=range(20))
            top10_acc = top_k_accuracy_score(true_labels, pred_probs_valid, k=10, labels=range(20))
            top20_acc = top_k_accuracy_score(true_labels, pred_probs_valid, k=20, labels=range(20))
            
            result = {
                'protein_name': os.path.splitext(os.path.basename(pdb_file))[0],
                'pdb_file': pdb_file,
                'total_residues': len(true_labels),
                'accuracy': accuracy,
                'top1_accuracy': top1_acc,
                'top3_accuracy': top3_acc,
                'top5_accuracy': top5_acc,
                'top10_accuracy': top10_acc,
                'top20_accuracy': top20_acc,
                'ts50_accuracy': ts50_acc,
                't500_accuracy': t500_acc,
                'true_labels': true_labels,
                'pred_probs': pred_probs_valid
            }
            
            return result
            
        except Exception as e:
            print(f"Error evaluating {pdb_file}: {e}")
            return None
    
    def calculate_multiple_proteins(self, pdb_files, show_progress=True):
        """Calculate T500 and TS50 for multiple proteins"""
        print(f"🧬 Calculating T500 and TS50 metrics for {len(pdb_files)} proteins...")
        
        if show_progress:
            pdb_files = tqdm(pdb_files, desc="Processing proteins")
        
        for pdb_file in pdb_files:
            result = self.calculate_t500_ts50_single_protein(pdb_file)
            if result:
                self.results.append(result)
        
        print(f"✅ Completed evaluation of {len(self.results)} proteins")
        return self.results
    
    def calculate_overall_metrics(self):
        """Calculate overall T500 and TS50 metrics"""
        if not self.results:
            print("No results available. Run evaluation first.")
            return None
        
        # Aggregate metrics
        total_residues = sum(r['total_residues'] for r in self.results)
        
        # Weighted averages (by number of residues)
        weighted_t500 = sum(r['t500_accuracy'] * r['total_residues'] for r in self.results) / total_residues
        weighted_ts50 = sum(r['ts50_accuracy'] * r['total_residues'] for r in self.results) / total_residues
        weighted_top1 = sum(r['top1_accuracy'] * r['total_residues'] for r in self.results) / total_residues
        weighted_top3 = sum(r['top3_accuracy'] * r['total_residues'] for r in self.results) / total_residues
        weighted_top5 = sum(r['top5_accuracy'] * r['total_residues'] for r in self.results) / total_residues
        weighted_top10 = sum(r['top10_accuracy'] * r['total_residues'] for r in self.results) / total_residues
        weighted_top20 = sum(r['top20_accuracy'] * r['total_residues'] for r in self.results) / total_residues
        
        # Simple averages
        simple_t500 = np.mean([r['t500_accuracy'] for r in self.results])
        simple_ts50 = np.mean([r['ts50_accuracy'] for r in self.results])
        simple_top1 = np.mean([r['top1_accuracy'] for r in self.results])
        simple_top3 = np.mean([r['top3_accuracy'] for r in self.results])
        simple_top5 = np.mean([r['top5_accuracy'] for r in self.results])
        
        self.overall_metrics = {
            'total_proteins': len(self.results),
            'total_residues': total_residues,
            
            # Weighted averages (more accurate for overall performance)
            'weighted_t500_accuracy': weighted_t500,
            'weighted_ts50_accuracy': weighted_ts50,
            'weighted_top1_accuracy': weighted_top1,
            'weighted_top3_accuracy': weighted_top3,
            'weighted_top5_accuracy': weighted_top5,
            'weighted_top10_accuracy': weighted_top10,
            'weighted_top20_accuracy': weighted_top20,
            
            # Simple averages
            'simple_t500_accuracy': simple_t500,
            'simple_ts50_accuracy': simple_ts50,
            'simple_top1_accuracy': simple_top1,
            'simple_top3_accuracy': simple_top3,
            'simple_top5_accuracy': simple_top5,
            
            # Per-protein statistics
            'min_t500': min(r['t500_accuracy'] for r in self.results),
            'max_t500': max(r['t500_accuracy'] for r in self.results),
            'std_t500': np.std([r['t500_accuracy'] for r in self.results]),
            'min_ts50': min(r['ts50_accuracy'] for r in self.results),
            'max_ts50': max(r['ts50_accuracy'] for r in self.results),
            'std_ts50': np.std([r['ts50_accuracy'] for r in self.results])
        }
        
        return self.overall_metrics
    
    def print_results(self):
        """Print comprehensive T500 and TS50 results"""
        if not self.overall_metrics:
            print("No overall metrics available. Run calculate_overall_metrics() first.")
            return
        
        print("\n" + "="*80)
        print("🧬 GCNdesign T500 and TS50 METRICS")
        print("="*80)
        
        print(f"\n📊 DATASET SUMMARY:")
        print(f"   Total proteins evaluated: {self.overall_metrics['total_proteins']}")
        print(f"   Total residues analyzed: {self.overall_metrics['total_residues']:,}")
        
        print(f"\n🎯 OVERALL PERFORMANCE (Weighted by residue count):")
        print(f"   ✓ T500 Accuracy: {self.overall_metrics['weighted_t500_accuracy']:.4f} ({self.overall_metrics['weighted_t500_accuracy']*100:.2f}%)")
        print(f"   ✓ TS50 Accuracy: {self.overall_metrics['weighted_ts50_accuracy']:.4f} ({self.overall_metrics['weighted_ts50_accuracy']*100:.2f}%)")
        print(f"   ✓ Top-1 Accuracy: {self.overall_metrics['weighted_top1_accuracy']:.4f} ({self.overall_metrics['weighted_top1_accuracy']*100:.2f}%)")
        print(f"   ✓ Top-3 Accuracy: {self.overall_metrics['weighted_top3_accuracy']:.4f} ({self.overall_metrics['weighted_top3_accuracy']*100:.2f}%)")
        print(f"   ✓ Top-5 Accuracy: {self.overall_metrics['weighted_top5_accuracy']:.4f} ({self.overall_metrics['weighted_top5_accuracy']*100:.2f}%)")
        print(f"   ✓ Top-10 Accuracy: {self.overall_metrics['weighted_top10_accuracy']:.4f} ({self.overall_metrics['weighted_top10_accuracy']*100:.2f}%)")
        print(f"   ✓ Top-20 Accuracy: {self.overall_metrics['weighted_top20_accuracy']:.4f} ({self.overall_metrics['weighted_top20_accuracy']*100:.2f}%)")
        
        print(f"\n📈 PER-PROTEIN STATISTICS:")
        print(f"   T500 - Min: {self.overall_metrics['min_t500']:.4f}, Max: {self.overall_metrics['max_t500']:.4f}, Std: {self.overall_metrics['std_t500']:.4f}")
        print(f"   TS50 - Min: {self.overall_metrics['min_ts50']:.4f}, Max: {self.overall_metrics['max_ts50']:.4f}, Std: {self.overall_metrics['std_ts50']:.4f}")
        
        # Performance interpretation
        t500 = self.overall_metrics['weighted_t500_accuracy']
        ts50 = self.overall_metrics['weighted_ts50_accuracy']
        
        print(f"\n🏅 PERFORMANCE ASSESSMENT:")
        print(f"   T500 Performance: {t500*100:.1f}%")
        print(f"   TS50 Performance: {ts50*100:.1f}%")
        
        # Compare with literature values
        print(f"\n📊 COMPARISON WITH LITERATURE:")
        print(f"   GCNdesign (3 layers): T500 ~53.49%, TS50 ~47.06%")
        print(f"   GCNdesign (4 layers): T500 ~53.78%, TS50 ~47.42%")
        print(f"   GCNdesign (5 layers): T500 ~53.64%, TS50 ~47.61%")
        print(f"   GCNdesign (6 layers): T500 ~53.41%, TS50 ~47.61%")
        print(f"   DenseCPD: T500 ~53.24%, TS50 ~46.74%")
        print(f"   ProDCoNN: T500 ~52.82%, TS50 ~50.71%")
        print(f"   SPROF: T500 ~42.20%, TS50 ~40.25%")
        print(f"   SPIN2: T500 ~40.69%, TS50 ~39.16%")
        
        # Performance level assessment
        if t500 > 0.53:
            t500_level = "EXCELLENT"
        elif t500 > 0.50:
            t500_level = "GOOD"
        elif t500 > 0.45:
            t500_level = "MODERATE"
        else:
            t500_level = "NEEDS IMPROVEMENT"
            
        if ts50 > 0.47:
            ts50_level = "EXCELLENT"
        elif ts50 > 0.45:
            ts50_level = "GOOD"
        elif ts50 > 0.40:
            ts50_level = "MODERATE"
        else:
            ts50_level = "NEEDS IMPROVEMENT"
        
        print(f"\n🏆 PERFORMANCE LEVELS:")
        print(f"   T500 Level: {t500_level}")
        print(f"   TS50 Level: {ts50_level}")
    
    def save_results(self, save_dir="t500_ts50_results"):
        """Save T500 and TS50 results to files"""
        os.makedirs(save_dir, exist_ok=True)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # Save detailed results
        results_df = pd.DataFrame(self.results)
        results_df.to_csv(os.path.join(save_dir, f't500_ts50_detailed_{timestamp}.csv'), index=False)
        
        # Save overall metrics
        with open(os.path.join(save_dir, f't500_ts50_overall_{timestamp}.pkl'), 'wb') as f:
            import pickle
            pickle.dump(self.overall_metrics, f)
        
        # Create comprehensive report
        report = f"""
GCNdesign T500 and TS50 Metrics Report
======================================
Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}

DATASET SUMMARY:
- Total proteins evaluated: {self.overall_metrics['total_proteins']}
- Total residues analyzed: {self.overall_metrics['total_residues']:,}

OVERALL PERFORMANCE (Weighted by residue count):
- T500 Accuracy: {self.overall_metrics['weighted_t500_accuracy']:.4f} ({self.overall_metrics['weighted_t500_accuracy']*100:.2f}%)
- TS50 Accuracy: {self.overall_metrics['weighted_ts50_accuracy']:.4f} ({self.overall_metrics['weighted_ts50_accuracy']*100:.2f}%)
- Top-1 Accuracy: {self.overall_metrics['weighted_top1_accuracy']:.4f} ({self.overall_metrics['weighted_top1_accuracy']*100:.2f}%)
- Top-3 Accuracy: {self.overall_metrics['weighted_top3_accuracy']:.4f} ({self.overall_metrics['weighted_top3_accuracy']*100:.2f}%)
- Top-5 Accuracy: {self.overall_metrics['weighted_top5_accuracy']:.4f} ({self.overall_metrics['weighted_top5_accuracy']*100:.2f}%)
- Top-10 Accuracy: {self.overall_metrics['weighted_top10_accuracy']:.4f} ({self.overall_metrics['weighted_top10_accuracy']*100:.2f}%)
- Top-20 Accuracy: {self.overall_metrics['weighted_top20_accuracy']:.4f} ({self.overall_metrics['weighted_top20_accuracy']*100:.2f}%)

PER-PROTEIN STATISTICS:
- T500 - Min: {self.overall_metrics['min_t500']:.4f}, Max: {self.overall_metrics['max_t500']:.4f}, Std: {self.overall_metrics['std_t500']:.4f}
- TS50 - Min: {self.overall_metrics['min_ts50']:.4f}, Max: {self.overall_metrics['max_ts50']:.4f}, Std: {self.overall_metrics['std_ts50']:.4f}

LITERATURE COMPARISON:
- GCNdesign (3 layers): T500 ~53.49%, TS50 ~47.06%
- GCNdesign (4 layers): T500 ~53.78%, TS50 ~47.42%
- GCNdesign (5 layers): T500 ~53.64%, TS50 ~47.61%
- GCNdesign (6 layers): T500 ~53.41%, TS50 ~47.61%
- DenseCPD: T500 ~53.24%, TS50 ~46.74%
- ProDCoNN: T500 ~52.82%, TS50 ~50.71%
- SPROF: T500 ~42.20%, TS50 ~40.25%
- SPIN2: T500 ~40.69%, TS50 ~39.16%

GENERATED FILES:
- t500_ts50_detailed_{timestamp}.csv - Detailed per-protein results
- t500_ts50_overall_{timestamp}.pkl - Overall metrics data
- t500_ts50_report_{timestamp}.txt - This report
"""
        
        with open(os.path.join(save_dir, f't500_ts50_report_{timestamp}.txt'), 'w') as f:
            f.write(report)
        
        print(f"💾 Results saved to: {save_dir}/")
        print(f"   📊 t500_ts50_detailed_{timestamp}.csv")
        print(f"   📁 t500_ts50_overall_{timestamp}.pkl")
        print(f"   📄 t500_ts50_report_{timestamp}.txt")

def main():
    """Main function to run T500 and TS50 calculation"""
    print("🧬 ProtGCN T500 and TS50 Metrics Calculator")
    print("=" * 60)
    
    # Initialize calculator
    calculator = T500TS50Calculator(device='cpu')
    
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
    
    # Calculate T500 and TS50
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
    
    print("\n✅ T500 and TS50 analysis complete!")

if __name__ == "__main__":
    main()
