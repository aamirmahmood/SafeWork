# This script will fuse a given sequence with given imu data. Also will plot it with ground truth and comparsion

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas
import copy
import yaml
import logging
import time
import random
import pandas as pd

from rotations import angle_normalize, rpy_jacobian_axis_angle, skew_symmetric, Quaternion


np.set_printoptions(threshold=sys.maxsize)

class Fusion():
    def __init__(self, config:dict) -> None:
        
        self.config = config

        try:
            self.cam_to_imu_tr_mat = self._load_transformation_matrix()
            self.imu_acceleration, self.imu_angular_velocity, self.imu_timestamp = self._load_imu_data()
            self.ground_truth, self.vo_estimates, self.ground_truth_timestamp, self.vo_estimate_timestamp = self._load_poses()
            self.srs_data_gt, self.srs_data_est, self.len_ssr_srs_est, self.bs_count = self._load_srs_data()
            self.var_imu_f = float(self.config["variance"]["imu_acceleration"])
            self.var_imu_w = float(self.config["variance"]["imu_angular_velocity"])
            self.var_gps  = float(self.config["variance"]["ground_truth"])
            self.var_vo = float(self.config["variance"]["visual_odometry"])
        except:
            logging.error(f"Error during initialization of fusion", exc_info=True)
            exit()
        
    def _load_transformation_matrix(self):

        T_from_cam_to_imu = None
        
        with open(self.config["transformation_matrix"]["kitti_velodyne_to_cam"], 'r') as f:
            calib = f.readlines()

        R = np.array([float(x) for x in calib[1].strip().split(' ')[1:]]).reshape((3, 3))
        t = np.array([float(x) for x in calib[2].strip().split(' ')[1:]])[:, None]

        T_velo_ref0 = np.vstack((np.hstack((R, t)), np.array([0, 0, 0, 1])))

        with open(self.config["transformation_matrix"]["kitti_imu_to_velodyne"], 'r') as f:
            calib = f.readlines()

        R = np.array([float(x) for x in calib[1].strip().split(' ')[1:]]).reshape((3, 3))
        t = np.array([float(x) for x in calib[2].strip().split(' ')[1:]])[:, None]

        T_imu_velo = np.vstack((np.hstack((R, t)), np.array([0, 0, 0, 1])))

        T_from_imu_to_cam = T_imu_velo @ T_velo_ref0
        T_from_cam_to_imu = np.linalg.inv(T_from_imu_to_cam)

        return T_from_cam_to_imu

    def _load_poses(self):
        
        ground_truth, visual_odometry = [[], [], []], [[], [], []]

        gt_file = self.config["data"]["ground_truth"]
        vo_file = self.config["data"]["visual_odometry_estimate"]

        gt_poses = pandas.read_csv(gt_file, delimiter=' ', header = None)
        vo_poses = pandas.read_csv(vo_file, delimiter=' ', header = None)

        for i in range(len(gt_poses)):
            mat = np.array(gt_poses.iloc[i]).reshape((3,4))
            s_mat = mat[:3].flatten()
            ground_truth[0].append(float(s_mat[3]))
            ground_truth[1].append(float(s_mat[7]))
            ground_truth[2].append(float(s_mat[11]))

        for i in range(len(vo_poses)):
            mat = np.array(vo_poses.iloc[i]).reshape((3,4))
            s_mat = mat[:3].flatten()
            visual_odometry[0].append(float(s_mat[3]))
            visual_odometry[1].append(float(s_mat[7]))
            visual_odometry[2].append(float(s_mat[11]))

        vo_in_imu_coord = []
        for x,y,z in zip(visual_odometry[0],visual_odometry[1],visual_odometry[2]):
            VO = np.array([x, y, z, 1])
            transformed = self.cam_to_imu_tr_mat @ VO
            vo_in_imu_coord.append([transformed[0], transformed[1], transformed[2]])

        gt_in_imu_coord = []
        for x,y,z in zip(ground_truth[0],ground_truth[1],ground_truth[2]):
            GT = np.array([x, y, z, 1])
            transformed = self.cam_to_imu_tr_mat @ GT
            gt_in_imu_coord.append([transformed[0], transformed[1], transformed[2]])

        vo_in_imu_coord = np.array(vo_in_imu_coord)
        gt_in_imu_coord = np.array(gt_in_imu_coord)

        vo_timestamp = self.imu_timestamp
        gt_timestamp = self.imu_timestamp

        return gt_in_imu_coord, vo_in_imu_coord, gt_timestamp, vo_timestamp
    
    def _load_srs_data(self):
        
        df = pd.read_excel('data/simResultsSRS_new.xlsx')
        df.dropna(how='all', inplace=True)
        reference_data = df['Reference'].to_numpy()
        estimated_data = df['Estimated'].to_numpy()
        reference_xyz = reference_data.reshape(-1, 3)  # shape: (N, 3)
        estimated_xyz = estimated_data.reshape(-1, 3)  # shape: (N, 3)
        ref_x = reference_xyz[:, 0]
        ref_y = reference_xyz[:, 1]
        ref_z = reference_xyz[:, 2]

        est_x = estimated_xyz[:, 0]
        est_y = estimated_xyz[:, 1]
        est_z = estimated_xyz[:, 2]
        
        with open("data/srs_nlos_los_new.csv", "r") as f:
            data = f.readlines()
        
        bs_count = []
        data = data[1:]
        for dat in data:
            dat = dat.split(",")
            bs_count.append(float(dat[3]))
        
        if not len(bs_count) == len(ref_x):
            print("Length of base station count and reference data are not equal.")
            exit()  
            
            
        with open("data/seq_09_gps_with_site_count_voronoi.csv", "r") as f:
            data = f.readlines()
        data = data[1:]
        base_stations_count = []
        for dat in data:
            dat = dat.split(",")
            base_stations_count.append(int(dat[5]))
               
        return [ref_x, ref_y, ref_z], [est_x, est_y, est_z], len(ref_x), base_stations_count

    def _load_imu_data(self):
        
        file = self.config["data"]["inertial_estimates"]
        imu = np.load(file)
        #imu_pd = pandas.DataFrame(imu, columns=['header.stamp', 'linear_acceleration.x', 'linear_acceleration.y', 'linear_acceleration.z', 'angular_velocity.x', 'angular_velocity.y','angular_velocity.z'])
        imu_acceleration = []
        imu_angular_velocity = []
        imu_time_stamps = []
        for i in imu:
            imu_acceleration.append([i[1], i[2], i[3]])
            imu_angular_velocity.append([i[4], i[5], i[6]])
            imu_time_stamps.append(i[0])
        
        return np.array(imu_acceleration), np.array(imu_angular_velocity), np.array(imu_time_stamps)

    def _measurement_update(self, sensor_var, p_cov_check, y_k, p_check, v_check, h_jac):
        # 3.1 Compute Kalman Gain
        r_cov = np.eye(3)*sensor_var
        k_gain = p_cov_check @ h_jac.T @ np.linalg.inv((h_jac @ p_cov_check @ h_jac.T) + r_cov)

        # 3.2 Compute error state
        error_state = k_gain @ (y_k - p_check)

        # 3.3 Correct predicted state
        p_hat = p_check + error_state[0:3]
        v_hat = v_check + error_state[3:6]

        # 3.4 Compute corrected covariance
        p_cov_hat = (np.eye(9) - k_gain @ h_jac) @ p_cov_check

        return p_hat, v_hat, p_cov_hat

    def _utility(self):
        
        imu_drop_out = float(self.config["utility"]["imu_drop_out"])
        vo_drop_out = float(self.config["utility"]["visual_odometry_drop_out"])

        downsampled_imu_data = self.imu_acceleration
        downsampled_imu_time = self.imu_timestamp
        downsampled_vo_data = self.vo_estimates
        downsampled_vo_time = self.vo_estimate_timestamp

        no_to_drop_imu = int(self.imu_acceleration.shape[0] * (1-imu_drop_out))
        no_to_drop_vo = int(self.vo_estimates.shape[0] * (1-vo_drop_out))

        indices_to_drop_imu = sorted(random.sample(range(self.imu_acceleration.shape[0]), no_to_drop_imu))
        indices_to_drop_vo = sorted(random.sample(range(self.vo_estimates.shape[0]), no_to_drop_vo))

        downsampled_imu_data = self.imu_acceleration[indices_to_drop_imu]
        downsampled_imu_time = self.imu_timestamp[indices_to_drop_imu]
        downsampled_vo_data = self.vo_estimates[indices_to_drop_vo]
        downsampled_vo_time = self.vo_estimate_timestamp[indices_to_drop_vo]

        self.imu_acceleration = downsampled_imu_data
        self.imu_timestamp = downsampled_imu_time
        self.vo_estimates = downsampled_vo_data
        self.vo_estimate_timestamp = downsampled_vo_time

        if self.config["utility"]["use_ground_truth"]:
            gt_drop_out = float(self.config["utility"]["ground_truth_drop_out"])
            no_to_drop_gt = int(self.ground_truth.shape[0] * (1-gt_drop_out))
            downsampled_gt_data = self.ground_truth
            downsampled_gt_time = self.ground_truth_timestamp
            indices_to_drop_gt = sorted(random.sample(range(self.ground_truth.shape[0]), no_to_drop_gt))
            downsampled_gt_data = self.ground_truth[indices_to_drop_gt]
            downsampled_gt_time = self.ground_truth_timestamp[indices_to_drop_gt]

            return downsampled_gt_data, downsampled_gt_time

        return None, None

    def fuse(self):

        gt_data, gt_time = None, None
        if self.config["utility"]["state"] == True:
            gt_data, gt_time = self._utility()

        # Constants
        g = np.array([0, 0, -9.81])  # gravity
        l_jac = np.zeros([9, 6])
        l_jac[3:, :] = np.eye(6)  # motion model noise jacobian
        h_jac = np.zeros([3, 9])
        h_jac[:, :3] = np.eye(3)  # measurement model jacobian

        # Initial values
        p_est = np.zeros([self.ground_truth.shape[0], 3])  # position estimates
        v_est = np.zeros([self.ground_truth.shape[0], 3])  # velocity estimates
        p_cov = np.zeros([self.ground_truth.shape[0], 9, 9])  # covariance matrices at each timestep

        p_est[0] = self.ground_truth[0]
        v_est[0] = np.array([0, 0, 0])
        p_cov[0] = np.zeros(9)  # covariance of estimate

        try:
            for k in range(1, self.ground_truth.shape[0]):  # start at 1 b/c we have initial prediction from gt
                delta_t = self.imu_timestamp[k] - self.imu_timestamp[k - 1]

                # 1. Update state with IMU inputs
                f_ns = (self.imu_acceleration[k - 1]) + g # calculate sum of forces
                p_check = p_est[k - 1, :] + delta_t*v_est[k - 1, :] + 0.5*(delta_t**2)*f_ns
                v_check = v_est[k - 1, :] + delta_t*f_ns

                # 1.1 Linearize the motion model and compute Jacobians
                f_jac = np.eye(9) # motion model jacobian with respect to last state
                f_jac[0:3, 3:6] = np.eye(3)*delta_t
                f_jac[3:6, 6:9] = -skew_symmetric(self.imu_acceleration[k - 1])*delta_t

                # 2. Propagate uncertainty
                q_cov = np.zeros((6, 6)) # IMU noise covariance
                q_cov[0:3, 0:3] = delta_t**2 * np.eye(3)*self.var_imu_f
                q_cov[3:6, 3:6] = 0 #delta_t**2 * np.eye(3)*self.var_imu_w
                p_cov_check = f_jac @ p_cov[k - 1, :, :] @ f_jac.T + l_jac @ q_cov @ l_jac.T

                # 3. Check availability of VO estimates
                if self.imu_timestamp[k] in self.vo_estimate_timestamp[0:self.vo_estimates.shape[0]]:
                    vo_i = list(self.vo_estimate_timestamp).index(self.imu_timestamp[k])
                    p_check, v_check, p_cov_check = self._measurement_update(self.var_vo, p_cov_check, self.vo_estimates[vo_i], p_check, v_check, h_jac)

                # # 3. Check availability of GPS updates
                # if gt_data is not None:
                #     if self.imu_timestamp[k] in gt_time:
                #         gt_i = list(gt_time).index(self.imu_timestamp[k])
                        
                #         p_check, v_check, p_cov_check = self._measurement_update(self.var_vo, p_cov_check, gt_data[gt_i], p_check, v_check, h_jac)

                # 4. Check avalability of PRS updates
                try:
                    if float(self.bs_count[k]) > 2:
                        p_check, v_check, p_cov_check = self._measurement_update(self.var_gps, p_cov_check, [self.srs_data_est[1][k], self.srs_data_est[0][k], self.srs_data_est[2][k]] , p_check, v_check, h_jac)
                except Exception as e:
                    print(e)
                    pass
                
                # Update states (save)
                p_est[k, :] = p_check
                v_est[k, :] = v_check
                p_cov[k, :, :] = p_cov_check
                
                if k == self.len_ssr_srs_est:
                    break
                
        except:
            logging.error("Error", exc_info=True)

        if self.config["results"]["plot"]:
            plt.figure()
            plt.plot(self.vo_estimates.T[1], self.vo_estimates.T[0], label='Visual Odometry', c='b', zorder=3)
            plt.plot(self.ground_truth.T[1], self.ground_truth.T[0], "--", c='k', label='Ground Truth', zorder=2)
            
            prs_data_plot_only_los_x = []
            prs_data_plot_only_los_y = []
            for i in range(len(self.srs_data_est[0])):
                if float(self.bs_count[i]) > 2:
                    prs_data_plot_only_los_x.append(self.srs_data_est[0][i])
                    prs_data_plot_only_los_y.append(self.srs_data_est[1][i])
            plt.scatter(prs_data_plot_only_los_x, prs_data_plot_only_los_y, label='SRS-LOS', c='r', zorder=4, marker='o')
            
            #plt.plot(self.prs_data_est[0], self.prs_data_est[1], label='PRS', c='r', zorder=4)
            plt.plot(p_est[:,1], p_est[:,0], label='ES-EKF Estimation', c='g', zorder=5)
            plt.xlabel("Distance travelled in X direction. (m)")
            plt.ylabel("Distance travelled in Y direction. (m)")
            plt.legend()
            plt.grid()
            
            if self.config["results"]["save_plot"]:
                plt.savefig(self.config["results"]["path_to_save"] + "/plot_data.png")
            else:
                plt.show()

            try:
                path = self.config["results"]["path_to_save"]
                if os.path.isdir(path):
                    ...
                else:
                    os.mkdir(path)
                np.save(path + "/fusion", p_est)
                np.save(path + "/vo_used", self.vo_estimates)
                if gt_data is not None:
                    np.save(path + "/gt_used", gt_data)
                else:
                    np.save(path + "/gt_used", self.ground_truth)
                
            except:
                logging.error("Error during saving data. ", exc_info=True)
                

            # Assuming self.p_est and self.ground_truth are NumPy arrays of shape (N, 2)
            # where each row represents a [y, x] coordinate (or [x, y] depending on your convention).

            # Compute the Euclidean error at each point
            error = np.linalg.norm(p_est - self.ground_truth, axis=1)
            # Only use the x and y coordinates for error calculation
            error = np.linalg.norm(p_est[:, :2] - self.ground_truth[:, :2], axis=1)

            # Create a new figure for the heat map
            plt.figure(figsize=(10, 8))

            # Use a scatter plot with the error as the color value
            sc = plt.scatter(p_est[:, 1], p_est[:, 0], c=error, cmap='hot', 
                            label='ES-EKF Estimation Error', zorder=5)

            # Optionally, overlay the other trajectories for reference
            plt.plot(self.vo_estimates.T[1], self.vo_estimates.T[0], label='Visual Odometry', c='b', zorder=3)
            plt.plot(self.ground_truth.T[1], self.ground_truth.T[0], "--", c='k', label='Ground Truth', zorder=2)
            plt.plot(self.srs_data_est[0], self.srs_data_est[1], label='PRS', c='r', zorder=4)

            plt.xlabel("Distance travelled in X direction (m)")
            plt.ylabel("Distance travelled in Y direction (m)")
            plt.title("Heat Map of ES-EKF Estimation Error")
            plt.colorbar(sc, label='Error (m)')
            plt.legend()
            plt.grid(True)

            if self.config["results"]["save_plot"]:
                plt.savefig(self.config["results"]["path_to_save"] + "/plot_heatmap.png")
            else:
                plt.show()
            
            
                with open(self.config["results"]["path_to_save"] + "/all_data.txt", "w") as f:
                    f.write("gt_x,gt_y,vo_x,vo_y,srs_x,srs_y,base_station,esk_x,esk_y\n")
                    esk_est_x = p_est[:, 1]
                    esk_est_y = p_est[:, 0]
                    for i in range(len(self.ground_truth.T[0])):
                        try:
                            f.write(f"{self.ground_truth.T[1][i]},{self.ground_truth.T[0][i]},{self.vo_estimates.T[1][i]},{self.vo_estimates.T[0][i]},{self.srs_data_est[0][i]},{self.srs_data_est[1][i]},{self.bs_count[i]},{esk_est_x[i]},{esk_est_y[i]}\n")
                        except:
                            f.write(f"{self.ground_truth.T[1][i]},{self.ground_truth.T[0][i]},{self.vo_estimates.T[1][-1]},{self.vo_estimates.T[0][-1]},{self.srs_data_est[0][i]},{self.srs_data_est[1][i]},{self.bs_count[i]},{esk_est_x[i]},{esk_est_y[i]}\n")

config_file = "config_sensor_fusion.yaml"
config = None
with open(config_file, "r") as f:
    config = yaml.safe_load(f)
    f.close()

logging.basicConfig(format='%(asctime)s - %(message)s', datefmt='%H:%M:%S', level=logging.INFO)

fusion = Fusion(config=config)
fusion.fuse()