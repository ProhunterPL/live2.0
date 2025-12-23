#!/bin/bash
# AWS Batch Setup Script
# Uruchom lokalnie (z AWS CLI skonfigurowanym)

set -e

REGION=${AWS_REGION:-us-east-1}
ACCOUNT_ID=$(aws sts get-caller-identity --query Account --output text)

echo "🚀 Live 2.0 - AWS Batch Setup"
echo "=============================="
echo "Region: $REGION"
echo "Account ID: $ACCOUNT_ID"
echo ""

# 1. S3 Bucket
echo "📦 Creating S3 bucket..."
aws s3 mb s3://live2-artifacts --region $REGION 2>/dev/null || echo "Bucket already exists"

# Lifecycle policy
cat > /tmp/lifecycle.json <<EOF
{
  "Rules": [
    {
      "Id": "DeleteOldArtifacts",
      "Status": "Enabled",
      "Expiration": {
        "Days": 90
      }
    }
  ]
}
EOF

aws s3api put-bucket-lifecycle-configuration \
  --bucket live2-artifacts \
  --lifecycle-configuration file:///tmp/lifecycle.json

echo "✅ S3 bucket created"

# 2. ECR Repository
echo "📦 Creating ECR repository..."
aws ecr create-repository \
  --repository-name live2-simulation \
  --region $REGION \
  2>/dev/null || echo "Repository already exists"

ECR_URI="$ACCOUNT_ID.dkr.ecr.$REGION.amazonaws.com/live2-simulation"
echo "✅ ECR repository: $ECR_URI"

# 3. IAM Role for Job Execution
echo "🔐 Creating IAM role for job execution..."
cat > /tmp/trust-policy.json <<EOF
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Effect": "Allow",
      "Principal": {
        "Service": "ecs-tasks.amazonaws.com"
      },
      "Action": "sts:AssumeRole"
    }
  ]
}
EOF

aws iam create-role \
  --role-name Live2JobRole \
  --assume-role-policy-document file:///tmp/trust-policy.json \
  2>/dev/null || echo "Role already exists"

# Attach S3 policy
cat > /tmp/s3-policy.json <<EOF
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Effect": "Allow",
      "Action": [
        "s3:PutObject",
        "s3:GetObject"
      ],
      "Resource": "arn:aws:s3:::live2-artifacts/*"
    },
    {
      "Effect": "Allow",
      "Action": [
        "logs:CreateLogStream",
        "logs:PutLogEvents"
      ],
      "Resource": "arn:aws:logs:*:*:log-group:/aws/batch/job"
    }
  ]
}
EOF

aws iam put-role-policy \
  --role-name Live2JobRole \
  --policy-name S3ArtifactsAccess \
  --policy-document file:///tmp/s3-policy.json

echo "✅ IAM role created: Live2JobRole"

# 4. IAM User for DO Orchestrator
echo "🔐 Creating IAM user for DO orchestrator..."
aws iam create-user --user-name live2-do-orchestrator 2>/dev/null || echo "User already exists"

cat > /tmp/batch-policy.json <<EOF
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Effect": "Allow",
      "Action": [
        "batch:SubmitJob",
        "batch:DescribeJobs",
        "batch:CancelJob"
      ],
      "Resource": "*"
    }
  ]
}
EOF

aws iam put-user-policy \
  --user-name live2-do-orchestrator \
  --policy-name BatchAccess \
  --policy-document file:///tmp/batch-policy.json

echo "✅ IAM user created: live2-do-orchestrator"
echo "⚠️  Create access key: aws iam create-access-key --user-name live2-do-orchestrator"

# 5. Batch Service Roles (if not exist)
echo "🔐 Checking Batch service roles..."

# AWSBatchServiceRole
if ! aws iam get-role --role-name AWSBatchServiceRole &>/dev/null; then
    echo "Creating AWSBatchServiceRole..."
    # This requires manual setup or use AWS managed policy
    echo "⚠️  Please create AWSBatchServiceRole manually in IAM Console"
fi

# aws-ec2-spot-fleet-role
if ! aws iam get-role --role-name aws-ec2-spot-fleet-role &>/dev/null; then
    echo "Creating aws-ec2-spot-fleet-role..."
    # This requires manual setup
    echo "⚠️  Please create aws-ec2-spot-fleet-role manually in IAM Console"
fi

echo ""
echo "✅ AWS infrastructure setup complete!"
echo ""
echo "📋 Next steps:"
echo "1. Create access key for live2-do-orchestrator:"
echo "   aws iam create-access-key --user-name live2-do-orchestrator"
echo "2. Setup Batch compute environment and job queue (see docs)"
echo "3. Build and push Docker image to ECR: $ECR_URI"
echo ""

