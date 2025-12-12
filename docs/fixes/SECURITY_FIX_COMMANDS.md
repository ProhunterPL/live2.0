# Security Fix - Remove SMTP Credentials from Git History

**CRITICAL:** Real SMTP password was exposed in commit `7aeb45c`.

---

## ✅ Immediate Actions (Already Done)

- ✅ Removed real password from all code files
- ✅ Replaced with placeholder `your-password#with#hash`
- ✅ Created security incident documentation

---

## 🔴 REQUIRED: Credential Rotation

**YOU MUST CHANGE THE SMTP PASSWORD IMMEDIATELY:**

1. Go to your email provider (Gmail/other)
2. Generate new App Password
3. Update `.env` file with new password
4. Test email alerts

**The exposed password `ADkEa32XLG#J` is compromised and must not be used.**

---

## 🧹 Git History Cleanup

Since commit `7aeb45c` was already pushed to GitHub, you need to clean git history.

### Option 1: BFG Repo-Cleaner (Recommended)

```bash
# 1. Download BFG: https://rtyley.github.io/bfg-repo-cleaner/
# 2. Create passwords.txt file:
echo "ADkEa32XLG#J==>your-password#with#hash" > passwords.txt

# 3. Clone a fresh copy of your repo
cd ..
git clone --mirror https://github.com/ProhunterPL/live2.0.git live2.0-clean.git

# 4. Run BFG
java -jar bfg.jar --replace-text passwords.txt live2.0-clean.git

# 5. Clean up
cd live2.0-clean.git
git reflog expire --expire=now --all
git gc --prune=now --aggressive

# 6. Force push (WARNING: Rewrites history!)
git push --force

# 7. Update your local repo
cd ../live2.0
git fetch origin
git reset --hard origin/main
```

### Option 2: git filter-branch

```bash
# Remove file from history
git filter-branch --force --index-filter \
  "git rm --cached --ignore-unmatch docs/plans/ENV_VARIABLES_GUIDE.md" \
  --prune-empty --tag-name-filter cat -- --all

# Clean up
git reflog expire --expire=now --all
git gc --prune=now --aggressive

# Force push
git push --force --all
```

### Option 3: GitHub Support (If you can't force push)

Contact GitHub support to remove the commit from public history.

---

## ⚠️ WARNINGS

1. **Force push rewrites history** - coordinate with team
2. **All collaborators** need to re-clone repo after cleanup
3. **Backup** your repo before cleanup
4. **The password is compromised** - must be changed regardless

---

## 📋 Checklist

- [x] Remove password from code
- [x] Remove password from documentation  
- [ ] **CHANGE SMTP PASSWORD** (CRITICAL)
- [ ] Update `.env` with new password
- [ ] Test email alerts
- [ ] Clean git history (if commit was pushed)
- [ ] Notify team about history rewrite
- [ ] Add pre-commit hooks to prevent future incidents

---

**Priority: CRITICAL - Change password immediately!**
