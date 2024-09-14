手动设置口腔打孔TCP（暂不用）

```c++
void ZYXtest::setOralTCP()
{
	dX = m_Controls.lineEdit_Oral_x->text().toDouble();
	dY = m_Controls.lineEdit_Oral_y->text().toDouble();
	dZ = m_Controls.lineEdit_Oral_z->text().toDouble();
	dRx = m_Controls.lineEdit_Oral_rx->text().toDouble();
	dRy = m_Controls.lineEdit_Oral_ry->text().toDouble();
	dRz = m_Controls.lineEdit_Oral_rz->text().toDouble();
	PrintTCP_qt(m_Controls.textBrowser, dX, dY, dZ, dRx, dRy, dRz);
	int nRet1 = HRIF_SetTCP(0, 0, dX, dY, dZ, dRx, dRy, dRz);
	//int nRet1 = HRIF_SetTCP(0, 0, dTcp_X, dTcp_Y, dTcp_Z, dTcp_Rx, dTcp_Ry, dTcp_Rz);
	PrintResult(m_Controls.textBrowser, nRet1, "TCP Set");
}
```

器械标定（来源于骁扬）需进一步debug

```c++
void DentalAccuracy::on_pushButton_calibrateDrill_Robot_clicked()
{
	// This function is aimed to calculate m_T_handpieceRFtoDrill

	// Step 1: calculate T_calibratorRFtoDrill from the
	// PointSet probe_head_tail_mandible or probe_head_tail_maxilla in the dataStorage
	if (GetDataStorage()->GetNamedNode("probe_head_tail_mandible") == nullptr ||
		GetDataStorage()->GetNamedNode("probe_head_tail_maxilla") == nullptr)
	{
		m_Controls.textBrowser->append("probe_head_tail_mandible or probe_head_tail_maxilla is missing!");
		return;
	}

	auto probe_head_tail_mandible = dynamic_cast<mitk::PointSet*>(GetDataStorage()->GetNamedNode("probe_head_tail_mandible")->GetData());
	auto probe_head_tail_maxilla = dynamic_cast<mitk::PointSet*>(GetDataStorage()->GetNamedNode("probe_head_tail_maxilla")->GetData());

	if (probe_head_tail_mandible->GetSize() != 2 || probe_head_tail_maxilla->GetSize() != 2)
	{
		m_Controls.textBrowser->append("probe_head_tail_mandible or probe_head_tail_maxilla is problematic!");
		return;
	}

	auto probe_head_tail = mitk::PointSet::New();

	if (m_Controls.radioButton_maxilla->isChecked())
	{
		probe_head_tail = probe_head_tail_maxilla;
	}
	else
	{
		probe_head_tail = probe_head_tail_mandible;
	}

	// 计算法兰到EndRF的转换矩阵T_FlangeRFToEndRF
	auto T_EndRFToCamera = vtkMatrix4x4::New();
	auto T_CameraToCalibrator = vtkMatrix4x4::New();
	auto T_FlangeToEndRF = vtkMatrix4x4::New();

	T_EndRFToCamera->DeepCopy(m_T_CameratoEndRF);
	T_EndRFToCamera->Invert();
	T_CameraToCalibrator->DeepCopy(m_T_cameraToCalibratorRF);
	T_FlangeToEndRF->DeepCopy(m_T_FlangeRFtoEndRF); 
	//auto tmpTrans = vtkTransform::New();
	//tmpTrans->Identity();
	//tmpTrans->PostMultiply();
	//tmpTrans->SetMatrix(T_CameraToEndRF);

	auto probe_head = probe_head_tail->GetPoint(0);
	auto probe_tail = probe_head_tail->GetPoint(1);

	Eigen::Vector3d z_probeInCalibratorRF;
	z_probeInCalibratorRF[0] = probe_tail[0] - probe_head[0];
	z_probeInCalibratorRF[1] = probe_tail[1] - probe_head[1];
	z_probeInCalibratorRF[2] = probe_tail[2] - probe_head[2];
	z_probeInCalibratorRF.normalize();

	Eigen::Vector3d z_std{ 0,0,1 };

	Eigen::Vector3d rotAxis = z_std.cross(z_probeInCalibratorRF);

	rotAxis.normalize();

	if (rotAxis.norm() < 0.00001) // in case the rotAxis becomes a zero vector
	{
		rotAxis[0] = 1;
		rotAxis[1] = 0;
		rotAxis[2] = 0;
	}

	double rotAngle = 180 * acos(z_std.dot(z_probeInCalibratorRF)) / 3.141592654;

	auto trans_calibratorRFtoDrill = vtkTransform::New();
	trans_calibratorRFtoDrill->Identity();
	trans_calibratorRFtoDrill->PostMultiply();
	trans_calibratorRFtoDrill->RotateWXYZ(rotAngle, rotAxis[0], rotAxis[1], rotAxis[2]);
	trans_calibratorRFtoDrill->Update();

	auto T_calibratorRFtoDrill = trans_calibratorRFtoDrill->GetMatrix();
	T_calibratorRFtoDrill->SetElement(0, 3, probe_head[0]);
	T_calibratorRFtoDrill->SetElement(1, 3, probe_head[1]);
	T_calibratorRFtoDrill->SetElement(2, 3, probe_head[2]);

	// for (int i{ 0 }; i < 16; i++)
	// {
	// 	m_T_calibratorRFtoDrill[i] = T_calibratorRFtoDrill->GetData()[i];
	// }

	memcpy_s(m_T_calibratorRFtoDrill, sizeof(double) * 16, T_calibratorRFtoDrill->GetData(), sizeof(double) * 16);


	m_Stat_calibratorRFtoDrill = true;

	// Step 2: Obtain the camera data and assemble the matrix:
	// T_handpieceRFtoDrill = (T_cameraTohandpieceRF)^-1 * T_cameraToCalibratorRF * T_calibratorRFtoDrill


	auto trans_FlangetoTCP = vtkTransform::New();
	trans_FlangetoTCP->Identity();
	trans_FlangetoTCP->PostMultiply();
	trans_FlangetoTCP->SetMatrix(T_calibratorRFtoDrill);
	trans_FlangetoTCP->Concatenate(T_CameraToCalibrator);
	trans_FlangetoTCP->Concatenate(T_EndRFToCamera);
	trans_FlangetoTCP->Concatenate(T_FlangeToEndRF);
	trans_FlangetoTCP->Update();

	// Todo: the matrix below should be averaged for a time span before being stored into m_T_handpieceRFtoDrill
	auto T_FlanegtoTCP = trans_FlangetoTCP->GetMatrix();

	memcpy_s(m_T_FlangeToTCP, sizeof(double) * 16, T_FlanegtoTCP->GetData(), sizeof(double) * 16);
	Eigen::Vector3d RX;
	Eigen::Vector3d RY;
	Eigen::Vector3d RZ;
	Eigen::Matrix3d Re;

	RX[0] = T_FlanegtoTCP->GetElement(0, 0);
	RX[1] = T_FlanegtoTCP->GetElement(0, 1);
	RX[2] = T_FlanegtoTCP->GetElement(0, 2);

	RY[0] = T_FlanegtoTCP->GetElement(1, 0);
	RY[1] = T_FlanegtoTCP->GetElement(1, 1);
	RY[2] = T_FlanegtoTCP->GetElement(1, 2);

	RZ[0] = T_FlanegtoTCP->GetElement(2, 0);
	RZ[1] = T_FlanegtoTCP->GetElement(2, 1);
	RZ[2] = T_FlanegtoTCP->GetElement(2, 2);

	Re << RX[0], RY[0], RZ[0],
		RX[1], RY[1], RZ[1],
		RX[2], RY[2], RZ[2];

	// 欧拉角顺序zyx
	Eigen::Vector3d eulerAngle = Re.eulerAngles(2, 1, 0);

	double tcp[6];
	tcp[0] = T_FlanegtoTCP->GetElement(3, 0);
	tcp[1] = T_FlanegtoTCP->GetElement(3, 1);
	tcp[2] = T_FlanegtoTCP->GetElement(3, 2);
	tcp[3] = vtkMath::DegreesFromRadians(eulerAngle(2));
	tcp[4] = vtkMath::DegreesFromRadians(eulerAngle(1));
	tcp[5] = vtkMath::DegreesFromRadians(eulerAngle(0));
	m_Controls.textBrowser->append("tcp x y z is:");
	m_Controls.textBrowser->append(QString::number(tcp[0]) + " /" + QString::number(tcp[1]) + "/ " + QString::number(tcp[2]));

	m_Controls.textBrowser->append("tcp rx ry rz is:");
	m_Controls.textBrowser->append(QString::number(tcp[3]) + " /" + QString::number(tcp[4]) + " /" + QString::number(tcp[5]));
	m_Stat_handpieceRFtoDrill = true;

	m_Controls.textBrowser->append("Robot Calibration Succeeded!");


}

```





开启手动拖动 直线打孔(此时拖动依据的坐标系将按照工具TCP执行)

```c++
void ZYXtest::openForceControl()
{
	//// 开启力控
	//int nState = 1;
	//int nRet = HRIF_SetForceControlState(0, 0, nState);
	//PrintResult(m_Controls.textBrowser, nRet, "Opening force control");
	bool bEnable = true;
	int nRet = HRIF_SetForceFreeDriveMode(0, 0, bEnable);
	selectDirectionForceControl();
	PrintResult(m_Controls.textBrowser, nRet, "Opening force control");
}
```

```c++
void ZYXtest::selectDirectionForceControl()
{
	//确定力控的坐标系
	int nMode = 1; //1:沿着TCP方向
	int nRet1 = HRIF_SetForceToolCoordinateMotion(0, 0, nMode);
	PrintResult(m_Controls.textBrowser, nRet1, "The force control coordinate system is set to TOOL");

	////// 定义力控策略
	//int nStrategy = 0;//0：定义为恒力模式  2：定义为越障模式
	//int nRet2 = HRIF_SetForceControlStrategy(0, 0, nStrategy);// 设置力控策略为恒力模式
	//PrintResult(m_Controls.textBrowser, nRet2, "Set the force control policy to constant force mode");

	// 定义力控自由驱动自由度状态（那个为1那个轴就是开启状态）
	int nX = 1; int nY = 1; int nZ = 1;
	int nRx = 1; int nRy = 0; int nRz = 0;
	int nRet3 = HRIF_SetFreeDriveMotionFreedom(0, 0, nX, nY, nZ, nRx, nRy, nRz);// 设置力控自由驱动自由度状态
	PrintResult(m_Controls.textBrowser, nRet3, "Set the direction of the force control is along the Z axis");

}
```

关闭手动拖动

```c++
void ZYXtest::closeForceControl()
{
	//关闭力控
	/*int nState = 0;
	int nRet = HRIF_SetForceControlState(0, 0, nState);
	PrintResult(m_Controls.textBrowser, nRet, "close force control");*/
	bool bEnable = false;
	int nRet = HRIF_SetForceFreeDriveMode(0, 0, bEnable);
	PrintResult(m_Controls.textBrowser, nRet, "close force control");
}

```

开启自由驱动（随意驱动）

```c++
void ZYXtest::OpenFreeDrag()
{
	int nRet = HRIF_GrpOpenFreeDriver(0, 0);
	PrintResult(m_Controls.textBrowser, nRet, "Open free drag");
}
```

关闭自由驱动

```c++
void ZYXtest::TurnOffFreeDrag()
{
	int nRet = HRIF_GrpCloseFreeDriver(0, 0);
	PrintResult(m_Controls.textBrowser, nRet, "Turn off free drag");
}
```

+++

### 机械臂导航定位

脊柱项目的导航执行

```c++
void ZYXtest::OnAutoPositionStart()
{
	// 确定目标线：会读取线数据导入到P2和P3中
	auto targetLinePoints = dynamic_cast<mitk::PointSet*>(m_Controls.mitkNodeSelectWidget_imageTargetLine_5->GetSelectedNode()->GetData());
    //将节点数据转换为mitk::PointSet类型，这个类型通常用于存储点集合
    
	auto targetPoint_0 = targetLinePoints->GetPoint(0); // TCP frame origin should move to this point
	auto targetPoint_1 = targetLinePoints->GetPoint(1);
       //从PointSet中获取两个点targetPoint_0和targetPoint_1，这些点可能定义了一条直线或平面，其中targetPoint_0是TCP的原点
	std::cout << "targetPoint_0: (" << targetPoint_0[0] << ", " << targetPoint_0[1] << ", " << targetPoint_0[2] << ")" << std::endl;
	std::cout << "targetPoint_1: (" << targetPoint_1[0] << ", " << targetPoint_1[1] << ", " << targetPoint_1[2] << ")" << std::endl;


	double targetPointUnderBase_0[3]{ 0 };
	double targetPointUnderBase_1[3]{ 0 };
       //用于存储变换后的目标点坐标

	for (int i{ 0 }; i < 20; i++)//手动滤波
	{
		//获取机械臂配准矩阵T_BaseToBaseRF
		vtkMatrix4x4* vtkT_BaseToBaseRF = vtkMatrix4x4::New();
		vtkT_BaseToBaseRF->DeepCopy(T_BaseToBaseRF);

		//获取T_BaseRFToCamera
		auto vtkT_BaseRFToCamera = vtkMatrix4x4::New();
		vtkT_BaseRFToCamera->DeepCopy(T_CamToBaseRF);
		vtkT_BaseRFToCamera->Invert();

		//获取T_CameraToPatientRF
		auto vtkT_CameraToPatientRF = vtkMatrix4x4::New();
		vtkT_CameraToPatientRF->DeepCopy(T_CamToPatientRF);

		//获取T_PatientRFToImage
		auto vtkT_PatientRFToImage = vtkMatrix4x4::New();
		vtkT_PatientRFToImage->DeepCopy(T_PatientRFtoImage_SPI);

		//计算T_BaseToImage
		vtkNew<vtkTransform> Transform;
		Transform->Identity();
		Transform->PostMultiply();
		Transform->SetMatrix(vtkT_PatientRFToImage);
		Transform->Concatenate(vtkT_CameraToPatientRF);
		Transform->Concatenate(vtkT_BaseRFToCamera);
		Transform->Concatenate(vtkT_BaseToBaseRF);
		Transform->Update();
		auto vtkT_BaseToImage = Transform->GetMatrix();
		PrintMatrix("T_BaseToImage", vtkT_BaseToImage->GetData());


//Tips: 为每个目标点创建一个矩阵TargetMatrix，并将其与T_BaseToImage组合，得到从基座到每个目标点的变换矩阵。
        
        
		//获取P0点的坐标
		auto TargetMatrix_0 = vtkMatrix4x4::New();
		TargetMatrix_0->SetElement(0, 3, targetPoint_0[0]);
		TargetMatrix_0->SetElement(1, 3, targetPoint_0[1]);
		TargetMatrix_0->SetElement(2, 3, targetPoint_0[2]);

		vtkNew<vtkTransform> Trans;
		Trans->Identity();
		Trans->PostMultiply();
		Trans->SetMatrix(TargetMatrix_0);
		Trans->Concatenate(vtkT_BaseToImage);
		Trans->Update();
		auto vtkT_BaseToTarget_0 = Trans->GetMatrix();

        
		//获取P1点的坐标
		auto TargetMatrix_1 = vtkMatrix4x4::New();
		TargetMatrix_1->SetElement(0, 3, targetPoint_1[0]);
		TargetMatrix_1->SetElement(1, 3, targetPoint_1[1]);
		TargetMatrix_1->SetElement(2, 3, targetPoint_1[2]);

		vtkNew<vtkTransform> Trans1;
		Trans1->Identity();
		Trans1->PostMultiply();
		Trans1->SetMatrix(TargetMatrix_1);
		Trans1->Concatenate(vtkT_BaseToImage);
		Trans1->Update();
		auto vtkT_BaseToTarget_1 = Trans1->GetMatrix();

		//计算20个点
		targetPointUnderBase_0[0] += vtkT_BaseToTarget_0->GetElement(0, 3);
		targetPointUnderBase_0[1] += vtkT_BaseToTarget_0->GetElement(1, 3);
		targetPointUnderBase_0[2] += vtkT_BaseToTarget_0->GetElement(2, 3);

		targetPointUnderBase_1[0] += vtkT_BaseToTarget_1->GetElement(0, 3);
		targetPointUnderBase_1[1] += vtkT_BaseToTarget_1->GetElement(1, 3);
		targetPointUnderBase_1[2] += vtkT_BaseToTarget_1->GetElement(2, 3);

	}
	//取平均
	targetPointUnderBase_0[0] = targetPointUnderBase_0[0] / 20;
	targetPointUnderBase_0[1] = targetPointUnderBase_0[1] / 20;
	targetPointUnderBase_0[2] = targetPointUnderBase_0[2] / 20;

	targetPointUnderBase_1[0] = targetPointUnderBase_1[0] / 20;
	targetPointUnderBase_1[1] = targetPointUnderBase_1[1] / 20;
	targetPointUnderBase_1[2] = targetPointUnderBase_1[2] / 20;
	//将累加的坐标除以循环次数（20），得到平均坐标，这可能用于减少误差或抖动(手动滤波)
    
	//获取机械臂的T_BaseToFlanger
	double dX1 = 0; double dY1 = 0; double dZ1 = 0;
	double dRx1 = 0; double dRy1 = 0; double dRz1 = 0;
	int nRet = HRIF_ReadActTcpPos(0, 0, dX1, dY1, dZ1, dRx1, dRy1, dRz1);

	auto tmpTrans = vtkTransform::New();
	tmpTrans->PostMultiply();
	tmpTrans->RotateX(dRx1);
	tmpTrans->RotateY(dRy1);
	tmpTrans->RotateZ(dRz1);
	tmpTrans->Translate(dX1, dY1, dZ1);
	tmpTrans->Update();
	vtkSmartPointer<vtkMatrix4x4> VTKT_BaseToFlanger = tmpTrans->GetMatrix();

	//借用法兰的X轴方向(第一列)
	Eigen::Vector3d currentXunderBase;
	currentXunderBase[0] = VTKT_BaseToFlanger->GetElement(0, 0);
	currentXunderBase[1] = VTKT_BaseToFlanger->GetElement(1, 0);
	currentXunderBase[2] = VTKT_BaseToFlanger->GetElement(2, 0);
	currentXunderBase.normalize();

	std::cout << "currentXunderBase: " << currentXunderBase << std::endl;

	//在Base坐标系下目标坐标系X轴的方向向量
	Eigen::Vector3d targetXunderBase;
	targetXunderBase[0] = targetPointUnderBase_1[0] - targetPointUnderBase_0[0];
	targetXunderBase[1] = targetPointUnderBase_1[1] - targetPointUnderBase_0[1];
	targetXunderBase[2] = targetPointUnderBase_1[2] - targetPointUnderBase_0[2];
	targetXunderBase.normalize();

    
	MITK_INFO << "targetXunderBase" << targetXunderBase<<std::endl;

	Eigen::Vector3d  rotationAxis;
	rotationAxis = currentXunderBase.cross(targetXunderBase);


	double rotationAngle;//定义旋转角
      //使用当前X轴和目标X轴计算旋转轴，然后根据两个向量的点积计算旋转角
	if (currentXunderBase.dot(targetXunderBase) > 0) //如果向量的内积大于0，cos大于0（锐角）
	{
		rotationAngle = 180 * asin(rotationAxis.norm()) / 3.1415926;//求向量的模长（sin[rotationAngle]）,再取反三角
	}
	else //如果向量的内积小于0，cos小于0（钝角）
	{
		rotationAngle = 180 - 180 * asin(rotationAxis.norm()) / 3.1415926;
	}

	vtkNew<vtkTransform> tmpTransform;
	tmpTransform->PostMultiply();
	tmpTransform->Identity();
	tmpTransform->SetMatrix(VTKT_BaseToFlanger);
	tmpTransform->RotateWXYZ(rotationAngle, rotationAxis[0], rotationAxis[1], rotationAxis[2]);//旋转角度，和旋转向量
	tmpTransform->Update();

	auto testMatrix = tmpTransform->GetMatrix();
	PrintMatrix("targetMatrix", testMatrix->GetData());

	Eigen::Matrix3d Re;

	Re << testMatrix->GetElement(0, 0), testMatrix->GetElement(0, 1), testMatrix->GetElement(0, 2),
		testMatrix->GetElement(1, 0), testMatrix->GetElement(1, 1), testMatrix->GetElement(1, 2),
		testMatrix->GetElement(2, 0), testMatrix->GetElement(2, 1), testMatrix->GetElement(2, 2);

	Eigen::Vector3d eulerAngle = Re.eulerAngles(2, 1, 0);
	double x = targetPointUnderBase_0[0];
	double y = targetPointUnderBase_0[1];
	double z = targetPointUnderBase_0[2];
	double rx = 180 * eulerAngle[2] / 3.1415;
	double ry = 180 * eulerAngle[1] / 3.1415;
	double rz = 180 * eulerAngle[0] / 3.1415;

	m_Controls.textBrowser->append("-------------------------------------------------------------------------------------------");
	m_Controls.textBrowser->append("TCP_move");
	m_Controls.textBrowser->append("dx=" + QString::number(x));
	m_Controls.textBrowser->append("dy=" + QString::number(y));
	m_Controls.textBrowser->append("dz=" + QString::number(z));
	m_Controls.textBrowser->append("dRx=" + QString::number(rx));
	m_Controls.textBrowser->append("dRy=" + QString::number(ry));
	m_Controls.textBrowser->append("dRz=" + QString::number(rz));
	/*("-------------------------------------------------------------------------------------------");*/

	dX = x;
	dY = y;
	dZ = z;
	dRx = rx;
	dRy = ry;
	dRz = rz;


	int nTargetPoint = HRIF_MoveL(0, 0, dX, dY, dZ, dRx, dRy, dRz, dJ1, dJ2, dJ3, dJ4, dJ5, dJ6, sTcpName, sUcsName,
		dVelocity, dAcc, dRadius, nIsSeek, nIOBit, nIOState, strCmdID);//机械臂移动

	if (nTargetPoint == 0) {
		m_Controls.textBrowser->append("MOVE OK");
	}
	else {
		m_Controls.textBrowser->append("MOVE FAILED");
	}
	m_Controls.textBrowser->append("-------------------------------------------------------------------------------------------");

}
```

DI机器人导航定位

```c++
bool ZYXtest::InterpretImageLine_Oral()
{
	//确定目标线：会读取线数据导入到P2和P3中
	auto targetLinePoints = dynamic_cast<mitk::PointSet*>(m_Controls.mitkNodeSelectWidget_imageTargetLine_6->GetSelectedNode()->GetData());
	auto targetPoint_0 = targetLinePoints->GetPoint(0); // TCP frame origin should move to this point
	auto targetPoint_1 = targetLinePoints->GetPoint(1);
	std::cout << "targetPoint_0: (" << targetPoint_0[0] << ", " << targetPoint_0[1] << ", " << targetPoint_0[2] << ")" << std::endl;
	std::cout << "targetPoint_1: (" << targetPoint_1[0] << ", " << targetPoint_1[1] << ", " << targetPoint_1[2] << ")" << std::endl;


	double targetPointUnderBase_0[3]{ 0 };
	double targetPointUnderBase_1[3]{ 0 };

	for (int i{ 0 }; i < 20; i++)
	{
		//获取机械臂配准矩阵T_BaseToBaseRF
		vtkMatrix4x4* vtkT_BaseToBaseRF = vtkMatrix4x4::New();
		vtkT_BaseToBaseRF->DeepCopy(T_BaseToBaseRF);


		//获取T_BaseRFToCamera
		auto vtkT_BaseRFToCamera = vtkMatrix4x4::New();
		vtkT_BaseRFToCamera->DeepCopy(T_CamToBaseRF);
		vtkT_BaseRFToCamera->Invert();

		//获取T_CameraToPatientRF
		auto vtkT_CameraToPatientRF = vtkMatrix4x4::New();
		vtkT_CameraToPatientRF->DeepCopy(T_CamToPatientRF);

		//获取T_PatientRFToImage
		auto vtkT_PatientRFToImage = vtkMatrix4x4::New();
		vtkT_PatientRFToImage->DeepCopy(T_PatientRFtoImage);

		//计算T_BaseToImage
		vtkNew<vtkTransform> Transform;
		Transform->Identity();
		Transform->PostMultiply();
		Transform->SetMatrix(vtkT_PatientRFToImage);
		Transform->Concatenate(vtkT_CameraToPatientRF);
		Transform->Concatenate(vtkT_BaseRFToCamera);
		Transform->Concatenate(vtkT_BaseToBaseRF);
		Transform->Update();
		auto vtkT_BaseToImage = Transform->GetMatrix();

		//获取P0点的坐标
		auto TargetMatrix_0 = vtkMatrix4x4::New();
		TargetMatrix_0->SetElement(0, 3, targetPoint_0[0]);
		TargetMatrix_0->SetElement(1, 3, targetPoint_0[1]);
		TargetMatrix_0->SetElement(2, 3, targetPoint_0[2]);

		vtkNew<vtkTransform> Trans;
		Trans->Identity();
		Trans->PostMultiply();
		Trans->SetMatrix(TargetMatrix_0);
		Trans->Concatenate(vtkT_BaseToImage);
		Trans->Update();
		auto vtkT_BaseToTarget_0 = Trans->GetMatrix();

		//获取P1点的坐标
		auto TargetMatrix_1 = vtkMatrix4x4::New();
		TargetMatrix_1->SetElement(0, 3, targetPoint_1[0]);
		TargetMatrix_1->SetElement(1, 3, targetPoint_1[1]);
		TargetMatrix_1->SetElement(2, 3, targetPoint_1[2]);

		vtkNew<vtkTransform> Trans1;
		Trans1->Identity();
		Trans1->PostMultiply();
		Trans1->SetMatrix(TargetMatrix_1);
		Trans1->Concatenate(vtkT_BaseToImage);
		Trans1->Update();
		auto vtkT_BaseToTarget_1 = Trans1->GetMatrix();

		//计算20个点
		targetPointUnderBase_0[0] += vtkT_BaseToTarget_0->GetElement(0, 3);
		targetPointUnderBase_0[1] += vtkT_BaseToTarget_0->GetElement(1, 3);
		targetPointUnderBase_0[2] += vtkT_BaseToTarget_0->GetElement(2, 3);

		targetPointUnderBase_1[0] += vtkT_BaseToTarget_1->GetElement(0, 3);
		targetPointUnderBase_1[1] += vtkT_BaseToTarget_1->GetElement(1, 3);
		targetPointUnderBase_1[2] += vtkT_BaseToTarget_1->GetElement(2, 3);

	}
	//取平均
	targetPointUnderBase_0[0] = targetPointUnderBase_0[0] / 20;
	targetPointUnderBase_0[1] = targetPointUnderBase_0[1] / 20;
	targetPointUnderBase_0[2] = targetPointUnderBase_0[2] / 20;

	targetPointUnderBase_1[0] = targetPointUnderBase_1[0] / 20;
	targetPointUnderBase_1[1] = targetPointUnderBase_1[1] / 20;
	targetPointUnderBase_1[2] = targetPointUnderBase_1[2] / 20;

	//获取机械臂的T_BaseToFlanger
	double dX1 = 0; double dY1 = 0; double dZ1 = 0;
	double dRx1 = 0; double dRy1 = 0; double dRz1 = 0;
	int nRet = HRIF_ReadActTcpPos(0, 0, dX1, dY1, dZ1, dRx1, dRy1, dRz1);

	auto tmpTrans = vtkTransform::New();
	tmpTrans->PostMultiply();
	tmpTrans->RotateX(dRx1);
	tmpTrans->RotateY(dRy1);
	tmpTrans->RotateZ(dRz1);
	tmpTrans->Translate(dX1, dY1, dZ1);
	tmpTrans->Update();
	vtkSmartPointer<vtkMatrix4x4> VTKT_BaseToFlanger = tmpTrans->GetMatrix();

	//借用法兰的X轴方向(第一列)
	Eigen::Vector3d currentXunderBase;
	currentXunderBase[0] = VTKT_BaseToFlanger->GetElement(0, 0);
	currentXunderBase[1] = VTKT_BaseToFlanger->GetElement(1, 0);
	currentXunderBase[2] = VTKT_BaseToFlanger->GetElement(2, 0);
	currentXunderBase.normalize();


	//std::cout << "currentxunderBase: " << currentXunderBase << std::endl;

	//在Base坐标系下目标坐标系X轴的方向向量
	Eigen::Vector3d targetXunderBase;
	targetXunderBase[0] = targetPointUnderBase_1[0] - targetPointUnderBase_0[0];
	targetXunderBase[1] = targetPointUnderBase_1[1] - targetPointUnderBase_0[1];
	targetXunderBase[2] = targetPointUnderBase_1[2] - targetPointUnderBase_0[2];
	targetXunderBase.normalize();

	//MITK_INFO << "targetXunderBase" << targetXunderBase;

	Eigen::Vector3d  rotationAxis;
	rotationAxis = currentXunderBase.cross(targetXunderBase);


	double rotationAngle;//定义旋转角
	if (currentXunderBase.dot(targetXunderBase) > 0) //如果向量的内积大于0，cos大于0（锐角）
	{
		rotationAngle = 180 * asin(rotationAxis.norm()) / 3.1415926;//求向量的模长（sin[rotationAngle]）,再取反三角
	}
	else //如果向量的内积小于0，cos小于0（钝角）
	{
		rotationAngle = 180 - 180 * asin(rotationAxis.norm()) / 3.1415926;
	}

	vtkNew<vtkTransform> tmpTransform;
	tmpTransform->PostMultiply();
	tmpTransform->Identity();
	tmpTransform->SetMatrix(VTKT_BaseToFlanger);
	tmpTransform->RotateWXYZ(rotationAngle, rotationAxis[0], rotationAxis[1], rotationAxis[2]);//旋转角度，和旋转向量
	tmpTransform->Update();

	auto testMatrix = tmpTransform->GetMatrix();

	Eigen::Matrix3d Re;

	Re << testMatrix->GetElement(0, 0), testMatrix->GetElement(0, 1), testMatrix->GetElement(0, 2),
		testMatrix->GetElement(1, 0), testMatrix->GetElement(1, 1), testMatrix->GetElement(1, 2),
		testMatrix->GetElement(2, 0), testMatrix->GetElement(2, 1), testMatrix->GetElement(2, 2);

	Eigen::Vector3d eulerAngle = Re.eulerAngles(2, 1, 0);
	double x = targetPointUnderBase_0[0];
	double y = targetPointUnderBase_0[1];
	double z = targetPointUnderBase_0[2];
	double rx = 180 * eulerAngle[2] / 3.1415;
	double ry = 180 * eulerAngle[1] / 3.1415;
	double rz = 180 * eulerAngle[0] / 3.1415;

	m_Controls.textBrowser->append("-------------------------------------------------------------------------------------------");
	m_Controls.textBrowser->append("TCP_move");
	m_Controls.textBrowser->append("dx=" + QString::number(x));
	m_Controls.textBrowser->append("dy=" + QString::number(y));
	m_Controls.textBrowser->append("dz=" + QString::number(z));
	m_Controls.textBrowser->append("dRx=" + QString::number(rx));
	m_Controls.textBrowser->append("dRy=" + QString::number(ry));
	m_Controls.textBrowser->append("dRz=" + QString::number(rz));
	("-------------------------------------------------------------------------------------------");

	dX = x;
	dY = y;
	dZ = z;
	dRx = rx;
	dRy = ry;
	dRz = rz;

	//机械臂移动
	int nTargetPoint = HRIF_MoveL(0, 0, dX, dY, dZ, dRx, dRy, dRz, dJ1, dJ2, dJ3, dJ4, dJ5, dJ6, sTcpName, sUcsName,
		dVelocity, dAcc, dRadius, nIsSeek, nIOBit, nIOState, strCmdID);

	PrintResult(m_Controls.textBrowser, nTargetPoint, "Move");

	return true;
}
```
