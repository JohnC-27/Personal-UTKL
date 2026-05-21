import ROOT
import array

# 1. Generate some random data
data_values = [ROOT.gRandom.Gaus(0, 1) for _ in range(100)]
data_array = array.array('d', data_values)

# 2. Initialize TKDE
# Arguments: n_events, data_array, xMin, xMax, options, rho (bandwidth scale)
kde = ROOT.TKDE(len(data_array), data_array, -5, 5, "", 1.0)

# 3. Plotting
canvas = ROOT.TCanvas("c1", "KDE Example", 800, 600)

# Extract as a TF1 for easy drawing/integration
kde_func = kde.GetFunction()
kde_func.SetLineColor(ROOT.kBlue)
kde_func.SetLineWidth(3)

# Draw a histogram underneath for comparison
hist = ROOT.TH1D("h1", "Data Distribution;X;Density", 30, -5, 5)
for v in data_values: hist.Fill(v)
hist.Scale(1.0 / hist.Integral("width")) # Normalize hist to compare with KDE

hist.Draw("HIST")
kde_func.Draw("SAME")

canvas.Update()
input()