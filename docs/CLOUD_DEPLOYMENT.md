#  RAPTOR v2.1.0 Cloud Deployment Guide

**Deploy RNA-seq Analysis to AWS, GCP, and Azure**

Run RAPTOR on cloud computing platforms for large-scale analyses, parallel processing, and cost-effective computing.

---

##  Table of Contents

1. [Overview](#overview)
2. [Why Use Cloud?](#why-use-cloud)
3. [Cost Estimation](#cost-estimation)
4. [AWS Deployment](#aws-deployment)
5. [Google Cloud Platform](#google-cloud-platform)
6. [Microsoft Azure](#microsoft-azure)
7. [Cost Optimization](#cost-optimization)
8. [Monitoring & Management](#monitoring--management)
9. [Troubleshooting](#troubleshooting)
10. [Best Practices](#best-practices)

---

##  Overview

RAPTOR v2.1.0 supports deployment to three major cloud platforms:

- **AWS (Amazon Web Services)** - Most mature ecosystem
- **GCP (Google Cloud Platform)** - Best pricing, great for ML
- **Azure (Microsoft Azure)** - Best for enterprise

### Cloud Features

✅ **Auto-scaling** - Spin up instances on demand  
✅ **Parallel processing** - Run multiple samples simultaneously  
✅ **Spot/preemptible instances** - 60-90% cost savings  
✅ **No maintenance** - No hardware to manage  
✅ **Pay per use** - Only pay for what you use  
✅ **Global access** - Run from anywhere  

---

##  Why Use Cloud?

### Use Cloud When:

✅ **Large projects** (50+ samples)  
✅ **No local compute** (no HPC cluster)  
✅ **Occasional analyses** (don't need dedicated hardware)  
✅ **Need fast turnaround** (scale up temporarily)  
✅ **Team collaboration** (share resources)  
✅ **Benchmarking** (compare all 8 pipelines)  

### Use Local/HPC When:

❌ **Small projects** (<20 samples)  
❌ **Have powerful hardware** (64GB+ RAM, 32+ cores)  
❌ **Frequent analyses** (daily use)  
❌ **Sensitive data** (cannot leave institution)  
❌ **Limited budget** (cloud costs add up)  
❌ **No internet** (poor connectivity)  

---

##  Cost Estimation

### Price Comparison (Approximate)

**Small Project (20 samples, 3 pipelines):**
| Provider | On-Demand | Spot/Preemptible | Time |
|----------|-----------|------------------|------|
| AWS | $15-25 | $5-10 | 4-6 hours |
| GCP | $12-20 | $4-8 | 4-6 hours |
| Azure | $18-28 | $6-12 | 4-6 hours |

**Medium Project (50 samples, ensemble):**
| Provider | On-Demand | Spot/Preemptible | Time |
|----------|-----------|------------------|------|
| AWS | $40-60 | $12-25 | 8-12 hours |
| GCP | $35-55 | $10-20 | 8-12 hours |
| Azure | $45-70 | $15-28 | 8-12 hours |

**Large Project (200 samples, benchmarking):**
| Provider | On-Demand | Spot/Preemptible | Time |
|----------|-----------|------------------|------|
| AWS | $150-250 | $50-100 | 24-48 hours |
| GCP | $130-220 | $45-90 | 24-48 hours |
| Azure | $180-280 | $60-120 | 24-48 hours |

** Tip:** Always use spot/preemptible instances for 60-90% savings!

### Cost Calculator

```bash
# Estimate costs for your project
raptor cloud estimate \
  --samples 50 \
  --pipelines 3 \
  --provider aws

# Output:
Estimated Cost: $15-25 (spot instances)
Estimated Time: 8-12 hours
Recommended Instance: r5.4xlarge
```

---

##  AWS Deployment

### Prerequisites

1. **AWS Account** - Sign up at https://aws.amazon.com
2. **AWS CLI** - Install command-line tools
3. **IAM Permissions** - S3, EC2, Batch access

### Quick Start

```bash
# 1. Configure AWS credentials
aws configure
# Enter: Access Key ID, Secret Access Key, Region

# 2. Deploy RAPTOR
raptor cloud deploy \
  --provider aws \
  --data s3://my-bucket/data/ \
  --output s3://my-bucket/results/ \
  --instance-type r5.4xlarge \
  --spot-instances

# 3. Monitor progress
raptor cloud status --job-id <job-id>

# 4. Download results
aws s3 sync s3://my-bucket/results/ ./local_results/
```

### Detailed Setup

#### Step 1: Install AWS CLI

```bash
# Linux/macOS
curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
unzip awscliv2.zip
sudo ./aws/install

# Or use pip
pip install awscli

# Verify installation
aws --version
```

#### Step 2: Configure Credentials

```bash
# Configure AWS CLI
aws configure

# You'll be prompted for:
AWS Access Key ID: AKIAIOSFODNN7EXAMPLE
AWS Secret Access Key: wJalrXUtnFEMI/K7MDENG/bPxRfiCYEXAMPLEKEY
Default region: us-east-1
Default output format: json

# Test configuration
aws s3 ls
```

#### Step 3: Create S3 Buckets

```bash
# Create bucket for data
aws s3 mb s3://my-raptor-data

# Create bucket for results
aws s3 mb s3://my-raptor-results

# Upload your data
aws s3 sync ./local_data/ s3://my-raptor-data/
```

#### Step 4: Set Up IAM Permissions

Create IAM policy (`raptor-policy.json`):

```json
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Effect": "Allow",
      "Action": [
        "s3:GetObject",
        "s3:PutObject",
        "s3:ListBucket"
      ],
      "Resource": [
        "arn:aws:s3:::my-raptor-data/*",
        "arn:aws:s3:::my-raptor-results/*"
      ]
    },
    {
      "Effect": "Allow",
      "Action": [
        "ec2:RunInstances",
        "ec2:TerminateInstances",
        "ec2:DescribeInstances",
        "ec2:CreateTags"
      ],
      "Resource": "*"
    },
    {
      "Effect": "Allow",
      "Action": [
        "batch:SubmitJob",
        "batch:DescribeJobs",
        "batch:TerminateJob"
      ],
      "Resource": "*"
    }
  ]
}
```

Apply policy:
```bash
aws iam create-policy \
  --policy-name RAPTORPolicy \
  --policy-document file://raptor-policy.json
```

#### Step 5: Launch Analysis

```bash
# Full command with all options
raptor cloud deploy \
  --provider aws \
  --region us-east-1 \
  --data s3://my-raptor-data/ \
  --output s3://my-raptor-results/ \
  --instance-type r5.4xlarge \
  --spot-instances \
  --max-spot-price 0.50 \
  --pipelines 1,3,4 \
  --threads 16 \
  --auto-shutdown \
  --notify your.email@example.com
```

### AWS Instance Types

**Memory-Optimized (Recommended):**
| Instance | vCPU | RAM | Price/hour | Best For |
|----------|------|-----|------------|----------|
| r5.large | 2 | 16 GB | $0.126 | Small projects |
| r5.xlarge | 4 | 32 GB | $0.252 | Testing |
| r5.2xlarge | 8 | 64 GB | $0.504 | Medium projects |
| r5.4xlarge | 16 | 128 GB | $1.008 | Large projects |
| r5.8xlarge | 32 | 256 GB | $2.016 | Very large |

**Compute-Optimized:**
| Instance | vCPU | RAM | Price/hour | Best For |
|----------|------|-----|------------|----------|
| c5.4xlarge | 16 | 32 GB | $0.68 | CPU-intensive |
| c5.9xlarge | 36 | 72 GB | $1.53 | Parallel |

** Use Spot Instances for 70% savings!**

### Monitoring AWS Jobs

```bash
# Check job status
raptor cloud status --provider aws --job-id <id>

# View logs
aws logs tail /aws/batch/raptor --follow

# Check costs so far
aws ce get-cost-and-usage \
  --time-period Start=2025-11-01,End=2025-11-20 \
  --granularity DAILY \
  --metrics BlendedCost
```

---

##  Google Cloud Platform

### Prerequisites

1. **GCP Account** - Sign up at https://cloud.google.com
2. **gcloud CLI** - Install command-line tools
3. **IAM Permissions** - Storage, Compute access

### Quick Start

```bash
# 1. Configure GCP
gcloud auth login
gcloud config set project YOUR_PROJECT_ID

# 2. Deploy RAPTOR
raptor cloud deploy \
  --provider gcp \
  --data gs://my-bucket/data/ \
  --output gs://my-bucket/results/ \
  --machine-type n2-highmem-16 \
  --preemptible

# 3. Monitor progress
raptor cloud status --provider gcp --job-id <job-id>

# 4. Download results
gsutil -m rsync -r gs://my-bucket/results/ ./local_results/
```

### Detailed Setup

#### Step 1: Install gcloud CLI

```bash
# Linux
curl https://sdk.cloud.google.com | bash
exec -l $SHELL

# macOS
brew install google-cloud-sdk

# Initialize
gcloud init

# Verify
gcloud --version
```

#### Step 2: Create Project

```bash
# Create new project
gcloud projects create raptor-analysis --name="RAPTOR RNA-seq"

# Set as default
gcloud config set project raptor-analysis

# Enable required APIs
gcloud services enable compute.googleapis.com
gcloud services enable storage-api.googleapis.com
```

#### Step 3: Create Storage Buckets

```bash
# Create bucket
gsutil mb -l us-central1 gs://my-raptor-data
gsutil mb -l us-central1 gs://my-raptor-results

# Upload data
gsutil -m rsync -r ./local_data/ gs://my-raptor-data/

# Set lifecycle policy (auto-delete old results)
cat > lifecycle.json << EOF
{
  "lifecycle": {
    "rule": [
      {
        "action": {"type": "Delete"},
        "condition": {"age": 90}
      }
    ]
  }
}
EOF

gsutil lifecycle set lifecycle.json gs://my-raptor-results
```

#### Step 4: Launch Analysis

```bash
raptor cloud deploy \
  --provider gcp \
  --project raptor-analysis \
  --zone us-central1-a \
  --data gs://my-raptor-data/ \
  --output gs://my-raptor-results/ \
  --machine-type n2-highmem-16 \
  --preemptible \
  --pipelines 1,3,4 \
  --auto-shutdown \
  --notify your.email@example.com
```

### GCP Machine Types

**High Memory (Recommended):**
| Machine Type | vCPU | RAM | Price/hour | Best For |
|--------------|------|-----|------------|----------|
| n2-highmem-4 | 4 | 32 GB | $0.24 | Small projects |
| n2-highmem-8 | 8 | 64 GB | $0.48 | Medium projects |
| n2-highmem-16 | 16 | 128 GB | $0.96 | Large projects |
| n2-highmem-32 | 32 | 256 GB | $1.92 | Very large |

**Preemptible (70% cheaper!):**
| Machine Type | vCPU | RAM | Price/hour | Notes |
|--------------|------|-----|------------|-------|
| n2-highmem-16 | 16 | 128 GB | $0.29 | Can be interrupted |

** Preemptible instances are perfect for RAPTOR!**

### Monitoring GCP Jobs

```bash
# List running instances
gcloud compute instances list

# Check instance details
gcloud compute instances describe INSTANCE_NAME

# View logs
gcloud logging read "resource.type=gce_instance"

# Check costs
gcloud beta billing accounts list
```

---

##  Microsoft Azure

### Prerequisites

1. **Azure Account** - Sign up at https://azure.microsoft.com
2. **Azure CLI** - Install command-line tools
3. **Resource Group** - Container for resources

### Quick Start

```bash
# 1. Configure Azure
az login
az account set --subscription "Your Subscription"

# 2. Deploy RAPTOR
raptor cloud deploy \
  --provider azure \
  --data https://mystorageaccount.blob.core.windows.net/data/ \
  --output https://mystorageaccount.blob.core.windows.net/results/ \
  --vm-size Standard_E16s_v3 \
  --spot-instance

# 3. Monitor progress
raptor cloud status --provider azure --job-id <job-id>

# 4. Download results
azcopy sync \
  "https://mystorageaccount.blob.core.windows.net/results/" \
  "./local_results/"
```

### Detailed Setup

#### Step 1: Install Azure CLI

```bash
# Linux
curl -sL https://aka.ms/InstallAzureCLIDeb | sudo bash

# macOS
brew install azure-cli

# Windows
# Download installer from https://aka.ms/installazurecliwindows

# Verify
az --version
```

#### Step 2: Login and Setup

```bash
# Login to Azure
az login

# List subscriptions
az account list --output table

# Set active subscription
az account set --subscription "Your Subscription Name"

# Create resource group
az group create \
  --name raptor-resources \
  --location eastus
```

#### Step 3: Create Storage Account

```bash
# Create storage account
az storage account create \
  --name raptorstorage \
  --resource-group raptor-resources \
  --location eastus \
  --sku Standard_LRS

# Get connection string
az storage account show-connection-string \
  --name raptorstorage \
  --resource-group raptor-resources

# Create containers
az storage container create \
  --name data \
  --account-name raptorstorage

az storage container create \
  --name results \
  --account-name raptorstorage

# Upload data
az storage blob upload-batch \
  --destination data \
  --source ./local_data/ \
  --account-name raptorstorage
```

#### Step 4: Launch Analysis

```bash
raptor cloud deploy \
  --provider azure \
  --resource-group raptor-resources \
  --location eastus \
  --data https://raptorstorage.blob.core.windows.net/data/ \
  --output https://raptorstorage.blob.core.windows.net/results/ \
  --vm-size Standard_E16s_v3 \
  --spot-instance \
  --pipelines 1,3,4 \
  --auto-shutdown \
  --notify your.email@example.com
```

### Azure VM Sizes

**Memory-Optimized:**
| VM Size | vCPU | RAM | Price/hour | Best For |
|---------|------|-----|------------|----------|
| Standard_E4s_v3 | 4 | 32 GB | $0.25 | Small projects |
| Standard_E8s_v3 | 8 | 64 GB | $0.50 | Medium projects |
| Standard_E16s_v3 | 16 | 128 GB | $1.01 | Large projects |
| Standard_E32s_v3 | 32 | 256 GB | $2.02 | Very large |

**Spot Instances (up to 90% savings!):**
| VM Size | vCPU | RAM | Price/hour | Max Price |
|---------|------|-----|------------|-----------|
| Standard_E16s_v3 | 16 | 128 GB | $0.10-0.30 | Variable |

### Monitoring Azure Jobs

```bash
# List VMs
az vm list --resource-group raptor-resources --output table

# Check VM status
az vm show \
  --resource-group raptor-resources \
  --name raptor-vm

# View activity log
az monitor activity-log list \
  --resource-group raptor-resources

# Check costs
az consumption usage list
```

---

##  Cost Optimization

### 8 Ways to Save Money

#### 1. Use Spot/Preemptible Instances (60-90% savings)

```bash
# AWS
raptor cloud deploy --spot-instances --max-spot-price 0.50

# GCP
raptor cloud deploy --preemptible

# Azure
raptor cloud deploy --spot-instance
```

#### 2. Auto-Shutdown When Complete

```bash
raptor cloud deploy --auto-shutdown
```

#### 3. Choose Right Region

**Cheaper regions:**
- AWS: us-east-1 (Virginia)
- GCP: us-central1 (Iowa)
- Azure: East US

**More expensive:**
- Any region with "2" (newer hardware)
- Asia-Pacific regions
- Europe (some)

#### 4. Use Appropriate Instance Size

```bash
# Don't over-provision!
# For 50 samples:
--instance-type r5.2xlarge  # Good ✅
--instance-type r5.8xlarge  # Overkill, wastes money ❌
```

#### 5. Compress Data

```bash
# Before upload
tar -czf data.tar.gz data/
# Saves storage costs and transfer time
```

#### 6. Delete Results After Download

```bash
# Set lifecycle policy to auto-delete after 30 days
# Or manually delete
aws s3 rm s3://my-raptor-results/ --recursive
```

#### 7. Use Reserved Instances (for frequent use)

If you run RAPTOR monthly:
- **AWS:** Reserved Instances (up to 72% savings)
- **GCP:** Committed Use Discounts (up to 57% savings)
- **Azure:** Reserved VM Instances (up to 72% savings)

#### 8. Monitor and Set Budgets

```bash
# AWS Budget Alert
aws budgets create-budget \
  --account-id 123456789 \
  --budget file://budget.json

# budget.json:
{
  "BudgetLimit": {
    "Amount": "100",
    "Unit": "USD"
  },
  "BudgetName": "RAPTOR Monthly Budget",
  "BudgetType": "COST",
  "TimeUnit": "MONTHLY"
}
```

---

##  Monitoring & Management

### Real-Time Monitoring

```bash
# Check job status
raptor cloud status --job-id <id>

# Watch progress
watch -n 30 'raptor cloud status --job-id <id>'

# Get detailed logs
raptor cloud logs --job-id <id> --follow
```

### Email Notifications

```bash
raptor cloud deploy \
  --notify your.email@example.com \
  --notify-on-complete \
  --notify-on-error
```

**You'll receive emails for:**
- ✅ Job started
- ✅ Job completed
- ❌ Job failed
- 💰 Cost threshold exceeded

### Dashboard Monitoring

```bash
# Start monitoring dashboard
raptor cloud dashboard \
  --provider aws \
  --port 8501

# View at http://localhost:8501
# Shows:
# - Running jobs
# - Cost tracking
# - Progress bars
# - Resource usage
```

---

##  Troubleshooting

### Issue: Permission Denied

**AWS:**
```bash
# Check IAM permissions
aws iam get-user-policy --user-name YOUR_USER --policy-name RAPTORPolicy

# Add missing permissions
aws iam attach-user-policy \
  --user-name YOUR_USER \
  --policy-arn arn:aws:iam::aws:policy/AmazonS3FullAccess
```

**GCP:**
```bash
# Check permissions
gcloud projects get-iam-policy raptor-analysis

# Add role
gcloud projects add-iam-policy-binding raptor-analysis \
  --member="user:you@example.com" \
  --role="roles/storage.admin"
```

**Azure:**
```bash
# Assign role
az role assignment create \
  --assignee you@example.com \
  --role "Storage Blob Data Contributor" \
  --scope "/subscriptions/SUBSCRIPTION_ID"
```

### Issue: Job Terminated (Spot Instance)

**Symptoms:**
```
Job terminated unexpectedly
Spot instance interrupted
```

**Solutions:**

1. **Enable checkpointing:**
```bash
raptor cloud deploy --checkpoint --resume-on-interrupt
```

2. **Use on-demand for critical jobs:**
```bash
raptor cloud deploy --on-demand  # No spot instances
```

3. **Set higher max price:**
```bash
raptor cloud deploy --spot-instances --max-spot-price 1.00
```

### Issue: Slow Data Transfer

**Solutions:**

1. **Use transfer acceleration (AWS):**
```bash
aws s3 cp data.tar.gz s3://bucket/ --acl bucket-owner-full-control \
  --storage-class INTELLIGENT_TIERING \
  --endpoint-url https://bucket.s3-accelerate.amazonaws.com
```

2. **Parallel transfer:**
```bash
# AWS
aws s3 sync data/ s3://bucket/ --parallel

# GCP
gsutil -m rsync -r data/ gs://bucket/

# Azure
azcopy sync data/ https://account.blob.core.windows.net/container/
```

3. **Compress first:**
```bash
tar -czf data.tar.gz data/
# Then upload compressed file
```

### Issue: High Costs

**Check current spend:**
```bash
# AWS
aws ce get-cost-and-usage \
  --time-period Start=2025-11-01,End=2025-11-20 \
  --granularity DAILY \
  --metrics BlendedCost

# GCP
gcloud beta billing accounts list

# Azure
az consumption usage list --subscription "Your Subscription"
```

**Stop everything:**
```bash
# AWS
raptor cloud terminate --provider aws --all

# GCP
gcloud compute instances delete --all

# Azure
az vm delete --resource-group raptor-resources --name raptor-vm
```

---

##  Best Practices

### Security

✅ **Use IAM roles** - Don't share access keys  
✅ **Encrypt data** - Enable encryption at rest  
✅ **Use VPC** - Isolate your instances  
✅ **Rotate credentials** - Regular key rotation  
✅ **Enable logging** - Track all access  
✅ **De-identify data** - Remove PHI before upload  

### Cost Management

✅ **Set budget alerts** - Get notified early  
✅ **Use spot instances** - 70% savings  
✅ **Auto-shutdown** - Don't leave running  
✅ **Delete old data** - Clean up regularly  
✅ **Right-size instances** - Don't over-provision  
✅ **Monitor usage** - Check costs daily  

### Data Management

✅ **Compress before upload** - Save bandwidth  
✅ **Use lifecycle policies** - Auto-delete old files  
✅ **Download results** - Don't rely on cloud storage  
✅ **Version your data** - Track what you uploaded  
✅ **Backup critical results** - Multiple locations  

### Performance

✅ **Choose nearby region** - Faster transfer  
✅ **Use SSD storage** - Faster I/O  
✅ **Parallel processing** - Multiple samples at once  
✅ **Monitor resources** - Ensure not bottlenecked  
✅ **Optimize instance type** - Balance cost vs speed  

---

##  Example Workflows

### Workflow 1: Small Project (Quick Test)

```bash
# 1. Upload data
aws s3 sync ./data/ s3://my-bucket/data/

# 2. Run single pipeline (fast)
raptor cloud deploy \
  --provider aws \
  --data s3://my-bucket/data/ \
  --output s3://my-bucket/results/ \
  --pipeline 3 \
  --instance-type r5.xlarge \
  --spot-instances \
  --auto-shutdown

# 3. Download results
aws s3 sync s3://my-bucket/results/ ./results/

# Cost: ~$5, Time: 2 hours
```

### Workflow 2: Medium Project (Recommended Workflow)

```bash
# 1. Upload data
gsutil -m rsync -r ./data/ gs://my-bucket/data/

# 2. Run ensemble (3 pipelines)
raptor cloud deploy \
  --provider gcp \
  --data gs://my-bucket/data/ \
  --output gs://my-bucket/results/ \
  --pipelines 1,3,4 \
  --ensemble \
  --machine-type n2-highmem-16 \
  --preemptible \
  --auto-shutdown \
  --notify me@example.com

# 3. Download results when notified
gsutil -m rsync -r gs://my-bucket/results/ ./results/

# Cost: ~$15, Time: 8 hours
```

### Workflow 3: Large Project (Full Benchmarking)

```bash
# 1. Upload data (compressed)
tar -czf data.tar.gz data/
aws s3 cp data.tar.gz s3://my-bucket/

# 2. Run all pipelines
raptor cloud deploy \
  --provider aws \
  --data s3://my-bucket/data.tar.gz \
  --output s3://my-bucket/results/ \
  --pipelines all \
  --benchmark \
  --instance-type r5.4xlarge \
  --spot-instances \
  --parallel 4 \
  --auto-shutdown \
  --notify me@example.com

# 3. Download comprehensive results
aws s3 sync s3://my-bucket/results/ ./results/

# Cost: ~$80, Time: 24 hours
```

---

## 🎉 Summary

Cloud deployment provides:
- ✅ **Scalability** - Handle any project size
- ✅ **Cost-effectiveness** - Pay only for what you use
- ✅ **Speed** - Parallel processing
- ✅ **Flexibility** - Choose optimal resources
- ✅ **Accessibility** - Run from anywhere
- ✅ **No maintenance** - Cloud provider handles hardware

**With spot/preemptible instances, cloud can be cheaper than maintaining local hardware!** 💰

---

##  Support

**Cloud deployment issues?**

1. Check [TROUBLESHOOTING.md](TROUBLESHOOTING.md)
2. Read [FAQ.md](FAQ.md) - Cloud section
3. GitHub Issues: https://github.com/AyehBlk/RAPTOR/issues
4. Email: ayehbolouki1988@gmail.com

**Platform-specific help:**
- AWS: https://docs.aws.amazon.com/
- GCP: https://cloud.google.com/docs
- Azure: https://docs.microsoft.com/azure/

---

**Author:** Ayeh Bolouki   
**Version:** 2.1.0  
**License:** MIT

---

*"Scale your RNA-seq analysis to the cloud!"* ☁️🚀
