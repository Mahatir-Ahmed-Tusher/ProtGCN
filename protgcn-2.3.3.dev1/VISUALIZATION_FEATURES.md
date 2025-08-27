# 🎨 ProtGCN Visualization Features

## Overview

ProtGCN now includes comprehensive visualization capabilities for both terminal and web interface predictions. The system generates professional-quality charts, graphs, and visual comparisons with the "ProtGCN" watermark.

## 🖼️ Generated Visualizations

### 1. **Sequence Comparison Chart**
- **Original vs Predicted sequences** side-by-side comparison
- **Color-coded accuracy**: Green for correct, red for incorrect predictions
- **Confidence scores** for each residue position
- **Amino acid properties** comparison (polar, nonpolar, acidic, basic)

### 2. **Accuracy & Performance Summary**
- **Overall accuracy pie chart**
- **Confidence score distribution histogram**
- **Amino acid frequency comparison** (original vs predicted)
- **Performance metrics summary** (accuracy, confidence levels)

### 3. **Confidence Heatmap**
- **Position-wise prediction probabilities** for all 20 amino acids
- **Color intensity** represents prediction confidence
- **Easy identification** of high/low confidence regions

## 🖥️ Terminal Usage

### Enhanced Prediction Script
```bash
# Basic prediction (text output only)
python scripts/gcndesign_predict_with_viz.py protein.pdb

# With visualizations (generates charts)
python scripts/gcndesign_predict_with_viz.py protein.pdb --visualize

# Custom output directory
python scripts/gcndesign_predict_with_viz.py protein.pdb --visualize --output-dir my_results

# With custom temperature
python scripts/gcndesign_predict_with_viz.py protein.pdb --visualize --temperature 0.5
```

### Generated Files
All visualization files are saved in the `predictions/` directory with timestamps:
- `protein_YYYYMMDD_HHMMSS_sequence_comparison.png`
- `protein_YYYYMMDD_HHMMSS_accuracy_summary.png`
- `protein_YYYYMMDD_HHMMSS_confidence_heatmap.png`

### Example Output
```bash
🧬 Starting ProtGCN prediction for 1ubq.pdb
🔧 Device: cpu, Temperature: 1.0

📊 PREDICTION RESULTS:
================================================================================
    1 M M:pred  0.703:M 0.047:Q 0.044:A 0.043:E 0.038:S
    2 Q T:pred  0.385:T 0.117:R 0.115:K 0.063:I 0.060:Q
    ...

🎨 Generating visualizations...
✅ Visualization complete!
📁 Files saved in: predictions/

📈 PREDICTION SUMMARY:
   Total residues: 76
   Accuracy: 51.32%
   Average confidence: 53.07%
   Predicted sequence: MTIYVADSDGTTYELEV...
   Original sequence:  MQIFVKTLTGKTITLEV...

🎯 Prediction completed!
📊 Check the 'predictions' folder for detailed visual analysis!
```

## 🌐 Web Interface Usage

### Automatic Visualization Generation
- **Automatic**: Visualizations are generated automatically for every prediction
- **No additional setup**: Just upload and predict as usual
- **Instant display**: Charts appear immediately after prediction
- **Download options**: Click download buttons to save charts

### Available Charts
1. **Sequence Comparison** (only when ground truth available)
2. **Accuracy & Performance Summary** (always available)
3. **Confidence Heatmap** (always available)

### Download Options
- **Individual downloads**: Click download button on each chart
- **High-quality PNG**: All downloads are 300 DPI publication-ready
- **Automatic naming**: Files named as `protgcn_[chart_type].png`

## 🎨 Visualization Features

### Color Schemes
- **Original sequence**: Blue (#2E86AB)
- **Predicted sequence**: Purple (#A23B72)
- **Correct predictions**: Green (#2E8B57)
- **Incorrect predictions**: Red (#DC143C)
- **Confidence scores**: Orange (#F18F01)

### Amino Acid Properties
- **Nonpolar**: Orange (A, V, L, I, M, F, W, P, G)
- **Polar**: Royal Blue (S, T, C, Y, N, Q)
- **Acidic**: Crimson (D, E)
- **Basic**: Forest Green (K, R, H)

### Watermarks
- **Position**: Bottom-right corner of all charts
- **Text**: "ProtGCN"
- **Style**: Semi-transparent, professional

## 📊 Chart Details

### Sequence Comparison Chart
**Components:**
- Top panel: Sequence accuracy (green/red bars)
- Middle panel: Confidence scores with thresholds
- Bottom panel: Amino acid property comparison

**Reading the chart:**
- Green bars = correct predictions
- Red bars = incorrect predictions
- Bar height in middle = confidence level
- Property colors show chemical similarity

### Accuracy Summary Chart
**Quadrant layout:**
1. **Top-left**: Accuracy pie chart
2. **Top-right**: Confidence distribution histogram
3. **Bottom-left**: Amino acid frequency comparison
4. **Bottom-right**: Performance metrics bars

### Confidence Heatmap
**Features:**
- Rows: 20 amino acids (A-Y)
- Columns: Protein positions
- Color intensity: Prediction probability
- Labels: Position/Original/Predicted

## 🔧 Technical Details

### Dependencies
```python
matplotlib>=3.10.5
seaborn>=0.13.2
numpy>=1.22.0
pandas>=2.3.2
```

### File Formats
- **Output**: PNG format
- **Resolution**: 300 DPI (publication quality)
- **Color space**: RGB
- **Background**: White

### Performance
- **Generation time**: 2-5 seconds per protein
- **File sizes**: 200KB - 2MB per chart
- **Memory usage**: ~50MB during generation

## 🚀 Advanced Usage

### Customizing Visualizations
The `ProtGCNVisualizer` class can be customized:

```python
from visualization import ProtGCNVisualizer

# Custom save directory
visualizer = ProtGCNVisualizer(save_dir="custom_output")

# Generate specific chart types
visualizer.create_sequence_comparison_chart(...)
visualizer.create_accuracy_summary_chart(...)
visualizer.create_confidence_heatmap(...)
```

### Batch Processing
```bash
# Process multiple proteins with visualizations
for protein in *.pdb; do
    python scripts/gcndesign_predict_with_viz.py "$protein" --visualize
done
```

### Integration with Other Tools
The visualization module can be imported and used independently:

```python
from visualization import ProtGCNVisualizer

# Create visualizer
viz = ProtGCNVisualizer()

# Generate web-ready base64 images
web_images = viz.generate_web_visualizations(results, summary, protein_name)

# Generate saved files
files = viz.generate_all_visualizations(results, summary, protein_name)
```

## 📝 Output Examples

### Terminal Output
```
📊 Generated 3 visualization files:
   ✓ predictions\1ubq_20250827_103818_sequence_comparison.png
   ✓ predictions\1ubq_20250827_103818_accuracy_summary.png
   ✓ predictions\1ubq_20250827_103818_confidence_heatmap.png
```

### Web Interface
- Interactive charts displayed inline
- Download buttons for each visualization
- Automatic generation with every prediction
- Mobile-responsive display

## 🎯 Benefits

### For Researchers
- **Publication-ready** figures
- **Comprehensive analysis** at a glance
- **Easy comparison** between predictions
- **Professional visualization** standards

### For Students
- **Visual learning** of protein sequence prediction
- **Clear confidence indicators**
- **Property-based understanding**
- **Interactive exploration**

### For Developers
- **Modular design** for easy integration
- **Customizable styling**
- **Multiple output formats**
- **Web and terminal compatibility**

---

**🧬 ProtGCN Visualization System - Making protein predictions visually compelling and scientifically rigorous! 🎨**
