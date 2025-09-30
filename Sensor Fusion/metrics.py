import numpy as np

def compute_ATE(gt, pred):
    
    ate = 0

    errors = []
    try:
        for idx,i in enumerate(pred):
        
            gt_xyz = gt[idx]
            pred_xyz = pred[idx]

            align_err = gt_xyz - pred_xyz
            errors.append(np.sqrt(np.sum(align_err ** 2)))
        
        ate = np.sqrt(np.mean(np.asarray(errors) ** 2))

    except IndexError:
        pass

    return ate

def compute_RPE(gt, pred, interval=1):
    errors = []
    for i in range(len(pred) - interval):
        delta_gt = np.linalg.norm(gt[i + interval] - gt[i])
        delta_pred = np.linalg.norm(pred[i + interval] - pred[i])
        errors.append(np.abs(delta_gt - delta_pred))
    return np.mean(errors)


# Example data loading and usage
# Ensure that the data arrays are correctly transposed if necessary
ground_truth = np.load("results_srs/gt_used.npy")
vo_predicted = np.load("results_srs/vo_used.npy")
kf_predicted = np.load("results_srs/fusion.npy")

# Pose alignment to first frame
kf_predicted_aligned = kf_predicted - kf_predicted[0]
vo_predicted_aligned = vo_predicted - vo_predicted[0]
ground_truth_aligned = ground_truth - ground_truth[0]

# Calculate ATE
ate_kf = compute_ATE(kf_predicted_aligned, ground_truth_aligned)
ate_vo = compute_ATE(kf_predicted_aligned, vo_predicted_aligned)
ate_vo_gt = compute_ATE(ground_truth_aligned, vo_predicted_aligned)
 
print(f"VO ATE: {ate_vo}")
print(f"KF ATE: {ate_kf}")
print(f"VOGR ATE: {ate_vo_gt}")


rpe_kf = compute_RPE(kf_predicted_aligned, ground_truth_aligned, 1)
rpe_vo = compute_RPE(kf_predicted_aligned, vo_predicted_aligned, 1)
rpe_vo_gt = compute_RPE(ground_truth_aligned, vo_predicted_aligned, 1)

print(f"VO RPE: {rpe_vo}")
print(f"KF RPE: {rpe_kf}")
print(f"VOGT RPE: {rpe_vo_gt}")